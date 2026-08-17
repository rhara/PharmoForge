# sequencealign 実装記録

このドキュメントは `sequencealign` 機能を再現するための仕様記録。

## 目的

ダウンロードした複数のPDB/CIF構造(結晶構造・AlphaFold予測構造)の蛋白配列を比較する。
「配列が本当に同一か」「どこに変異・構築上の違いがあるか」を素早く確認できるようにする。
ユーザーからの依頼: 「ダウンロードした蛋白の構造データのシーケンスを比較したい」。

## 検討・決定事項

- 実装範囲は「その場限りのスクリプト」ではなく、`pf`の正式なサブコマンド(`pf sequence-align`)として
  実装する方針をユーザーとの選択式確認で決定(他コマンドと同様、PROMPT記録・テスト・README付き)。
- コマンド名はユーザー指定により`sequence-align`(`pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a`
  のように使いたいという要望)。
- 配列アラインメントの実装方針についてユーザーに確認したところ「ProDyの方が好き」との回答。
  調査の結果、ProDy自体の配列アラインメント機能(`matchChains`等)も内部でBiopython
  (`Bio.Align.PairwiseAligner`)を呼ぶラッパーであることが判明。そのため:
  - 構造(Atomic)同士の比較は、ユーザーの意向を汲みProDyの`matchChains`をラップする
    `src/structcompare`を実装(呼び出しはProDy APIに閉じる)。
  - 後日追加した「構造を伴わない任意配列を基準にしたい」という要望(下記)に対しては、
    比較対象の一方に構造(Atomic)が存在しないためProDyでは対応できず、Biopythonの
    `PairwiseAligner`を直接使う`src/seqalign`を別途実装(この場合のみBiopython直接利用)。

## `pwalign`に関する重大な既知の不具合(実データで発見)

開発中、ProDy `matchChains`の`pwalign=True`(配列アラインメントへのフォールバック)について、
実データ(AlphaFold全長モデル vs 結晶構造のドメインのみの構築物)で重大な不具合を発見した:

- 返される`AtomMap`の残基番号(`getResnums()`)が実際の対応関係を反映しない
  (例: AF全長モデルの残基1が結晶構造の残基449に対応すると報告される等)。
- 検証: 重ね合わせ後のRMSDを計算すると数十Å(例: 38.6Å)に達し、対応が完全に破綻していることを確認。
- 一方、`pwalign=False`(残基番号+残基名による直接対応付け)は同じ入力で正しく機能する
  (重ね合わせ後RMSD 0.83Å)。返される%identity/%overlapの値自体は`pwalign`の有無によらず
  一致する(検証済み)。
- さらに`pwalign=True`は入力の組み合わせによっては内部の`getAlignedMatch`で`StopIteration`を
  送出して丸ごと失敗することも確認(例: `6LU7_ab` vs `8UPW_ab`)。

このため`structcompare`では:
- `match_chains()`(pairwise identity一覧)は`pwalign=True`を使う(情報提供目的、
  %identity/%overlapのみ使用。StopIterationはキャッチしてそのペアをスキップ)。
- `find_substitutions()`(残基単位の置換検出)は`pwalign=False`のみを使い、対応が
  取れない場合は`matched=False`を返す(既存の`structfit --method number`と同じ前提)。

またチェーンの一方が蛋白でない(残基0件)場合に自明な100%/100%マッチが返ることも確認し、
`_match()`内で`numAtoms() > 0`のものだけを残すようフィルタしている。

## 実データでの検証結果

- CDK2(P24941)のAlphaFoldモデルと結晶構造(1AQ1, 1HCL)は完全に配列一致(置換0)を確認。
- BRAF(P15056)のAlphaFoldモデルと結晶構造4MNFを比較すると、既知の発がん性変異**V600E**が
  正しく1件のみ検出された(実際の生物学的事実と一致する検証)。
- BRAF結晶構造3OG7では14箇所の相違を検出(結晶化用の構築物側の変異と見られる)。

## 欠損領域(ギャップ)の可視化

置換一覧のみでは、CDK2の1AQ1_ab/1HCL_aで実際に見えていた欠損領域(T-loop等の可動領域が
結晶構造で解析できていない)が分からないという指摘をユーザーから受け、追加実装した
(「変異だけでなく欠落しているギャップをみたい(PDBではギャップがかなりありますので)」)。

- `ラベル:チェーンID`基準(構造間比較): `structcompare.SubstitutionReport`に
  マッチしたチェーンID(`chain_id_b`)を追加し、`sequencealign/report.py`側で基準チェーンと
  対象チェーンの残基番号集合の差分(set演算)を取ってギャップとして表示する
  (`structcompare`自体は変更せず、`seqextract`で既に取得済みの`ChainSequence.resnums`を再利用)。
- アミノ酸配列基準(配列アラインメント): `seqalign.align_to_reference()`が返す型を
  `list[SequenceSubstitution]`から`AlignmentResult`(`identity`/`coverage`/`substitutions`/`gaps`)
  に拡張。アラインメントのブロック境界(`alignment.aligned`)から、基準配列側にのみ存在する
  領域(`kind="deletion"`)とquery側にのみ存在する領域(`kind="insertion"`、発現タグ等を想定)を
  検出する。
- 実データ検証(CDK2): `1AQ1_ab`で欠損36-43, 149-161、`1HCL_a`で欠損37-40を検出
  (いずれもCDK2の既知の可動領域(活性化ループ周辺)と符合)。
- 実データ検証(BRAF、`--reference P15056_AF:A`): `4MNF_ac`で欠損1-448, 601-615, 721-766
  (全長のうちキナーゼドメイン以外+活性化ループの一部が結晶構造で不可視)を検出。

## 整列表示(残基番号ベースの横並び表示)

ユーザーからの追加要望: 「シーケンスを1行100残基程度でタンパク質横断的にマッチしたシーケンスを
並べたい」。従来のFASTAセクションは各構造の配列を個別に(残基番号の対応を無視して)羅列するのみで、
構造間の対応関係が視覚的に分かりにくかった。

- `format_alignment_block()`(`src/sequencealign/report.py`)を追加。配列アラインメントは
  行わず、既存の前提(構造間でPDBの残基番号が揃っている)に従い、全構造・全チェーンの残基番号の
  和集合を軸として各配列を並べる(観測されていない残基は`-`で埋める)。`width`(既定100)残基
  ごとに`-- 開始-終了 --`の見出しで折り返すブロック形式(Clustal風)。
- レポートの新セクション`== 整列表示(残基番号ベース) ==`として常時出力(Pairwise identityと
  基準配列に対する置換の間に配置)。
- 異なる蛋白の構造を混在させた場合は意味のない結果になる点をREADME/docstringに明記
  (`structfit`/`align-view --method number`と同じ既存の前提を踏襲するのみで、新たな検証ロジックは
  追加していない)。

### ルーラー(位置番号目盛り)の追加と桁欠け不具合の修正

ユーザーからの追加要望: 「120 / | のように各ブロックごとに位置番号の目盛りを表示してほしい」。
`_format_ruler()`を追加し、各ブロックの直上に10残基ごとの位置番号(数字の行)と`|`(目盛り行)の
2行を表示するようにした。

実装当初、ブロック左端に近い目盛り(例: ブロックがresnum=449から始まる場合、最初の10の倍数である
450の目盛り)では、3桁の数字("450")のうち先頭の桁がブロック範囲外(負のインデックス)にあたる
ため描画時に捨てられ、**末尾の'0'だけが見える**(実際の値が誤って読める)不具合があった。
ユーザーから「各ブロックの一番左の値が常に0に見える」という指摘を受けて発見。

最初の修正では、数字全体がブロック内に収まらない目盛りをまるごと省略する方針にしたが、これだと
`--width`が100等キリの良い数字で、かつ構造の残基番号がキリの良い数字から始まらない場合
(例: 449始まり、`--width 100`だと2ブロック目が549始まり)、**2ブロック目以降のブロック先頭付近の
目盛りが毎回表示されない**という新たな指摘を受けた(「桁幅が10の倍数のとき2ブロック目から
ブロック先頭の桁が表示されません」)。これを受けて最終的な方針に変更: 数字が右揃えでは収まらない
目盛りは、`|`は本来の列のまま、数字だけをブロック左端(列0)に寄せて全体を表示する
(隣接する目盛り同士は10列以上離れており数字は最大4桁のため、この寄せによる重なりは生じない)。
これにより値が欠落することも誤読することもなくなる。回帰テストを`tests/sequencealign/test_report.py`の
`test_format_ruler_left_aligns_number_that_would_overflow_left_edge`に追加。

ルーラー導入後、ユーザーから「ブロック先頭の`-- nnn-nnn --`表示はもはやいらない」との指摘を受け、
`format_alignment_block()`から各ブロックの見出し行(`-- 開始-終了 --`)を削除した(ルーラー自体が
範囲を示すため冗長だった)。ブロックは空行のみで区切られる。関連テストのアサーションも
見出し文字列ではなく配列内容・ブロック数ベースに更新した。

### `--width`オプション

ユーザー要望: 「桁は100をデフォルトとして`--width 160`のように指定できると嬉しい」。
`format_alignment_block()`/`build_report()`に`width`/`align_width`引数を追加し、
`pf sequence-align --width <残基数>`(既定`DEFAULT_ALIGN_WIDTH=100`)で折り返し幅を指定できる
ようにした。

## コマンド出力の英語化

ユーザー要望: 「出力はすべて英語でお願いします」→ 範囲を確認したところ「pf sequence-alignの
レポート出力のみ」。続けて「--helpも含めてすべて英語にしてください」と追加指示があったため、
`src/sequencealign/report.py`が生成するレポート本文(セクション見出し・「no substitutions」
「gaps:」等)と`src/sequencealign/cli.py`の`--help`テキスト(docstring・各オプションのhelp文字列)を
すべて英語にした。README/PROMPT.md等のドキュメント、および他コマンド(`align-view`等)は
PharmoForgeの既存方針通り日本語のまま変更していない。

## `--reference`オプションの仕様

`--reference`は2通りの指定方法をサポートする(ユーザーからの追加要望により、後日拡張):

1. `ラベル:チェーンID`(例: `P24941_AF:A`): `--indir`等で読み込んだ構造の1チェーンを基準にする。
   残基番号ベースの対応付け(`structcompare.find_substitutions`)を用いる。
2. アミノ酸配列(1文字表記、コロンを含まない文字列): 構造を伴わない任意配列
   (UniProt正規配列やユーザーが直接貼り付けた配列)を基準にできるようにしたいという要望。
   例: `pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab 1HCL_a`。
   この場合は`seqalign.align_to_reference()`(Biopython `PairwiseAligner`によるグローバル
   アラインメント)を用いる。基準配列内の位置(`ref_pos`)と構造側の実際の残基番号
   (`query_resnum`)は一致するとは限らない(基準配列がその構造とUniProt番号で厳密に
   対応する保証がない)ため、出力では両方を明示する(例: `V600E(構造残基番号=600)`)。

判定方法: `reference`文字列に`:`が含まれるかどうかで機械的に判別する(構造ラベルは
`--indir`解決の都合上コロンを含み得ないため曖昧さはない)。コロンを含まない場合は
1文字アミノ酸コードの正規表現(標準20種+曖昧/非標準コード)で妥当性を検証し、
一致しなければエラーにする。

## CLI仕様

(`--reference`は後日廃止された。現在のCLI仕様は下記「FASTA入力対応と`--reference`オプションの
廃止」節、および[README.md](README.md)を参照。以下は導入当時の記録。)

```
pf sequence-align <構造ファイル1> [<構造ファイル2> ...] [--reference <基準>] [--width <残基数>] [--output <出力ファイル>]
```

- 構造ファイルの指定は`align-view`と同じトークン解決(`--indir`、拡張子省略時の自動補完)を
  共有する(下記「`--indir`解決ロジックの共通化」参照)。
- `--reference`省略時は整列表示に基準配列の行を加えない。
- `--width`(既定`DEFAULT_ALIGN_WIDTH=100`)は整列表示セクションの折り返し残基数。
- `--output`/`-o`省略時は標準出力にレポートを出す。

### 出力レポートの構成(`src/sequencealign/report.py`)

レポート本文はすべて英語(下記「コマンド出力の英語化」参照)。2セクションのみ
(下記「出力セクションの絞り込み」参照。当初はFASTA配列一覧・基準配列に対する置換一覧の
セクションもあったが、ユーザー要望により削除した)。

1. `== Pairwise identity ==`: 全構造の組み合わせについて`structcompare.match_chains()`の結果を一覧化。
2. `== Alignment (by residue number) ==`: 残基番号ベースの整列表示(上記「整列表示」参照)。
   `--reference`にアミノ酸配列を直接指定した場合は、`reference`行としてこのブロックにも
   加わる(上記「整列表示への基準配列の追加」参照)。`--reference`指定時は`_validate_reference()`
   により妥当性を検証する(存在しないラベル/チェーンや不正な配列文字列は`ValueError`。
   下記「出力セクションの絞り込み」参照)。

### 整列表示への基準配列の追加

上記対応の直後、ユーザーから「レファレンスのsequenceを、比較ブロックに出力してほしい。
`== Substitutions relative to reference ==`のセクションに入りません」との指摘を受けた。
`--reference`を直接アミノ酸配列で指定した場合、その配列は「Substitutions」セクションには
出るが「Alignment(比較ブロック)」セクションには出ておらず、他の構造と横並びで見比べられない
という指摘。

`format_alignment_block()`に`reference`引数を追加し、`reference`がコロンを含まない
(=アミノ酸配列直接指定)場合、基準配列の1文字目を残基番号1として扱った`reference`行を
整列表示のエントリ先頭に加えるようにした(基準配列がUniProt正規配列等、通常残基1から
始まる前提。既存の「構造間でPDBの残基番号が揃っている」前提の延長)。`reference`が
`ラベル:チェーンID`(構造)の場合は、対応するチェーンがすでに`structures`側の行として
含まれているため、重複追加はしない(`format_alignment_block(structures)`と
`format_alignment_block(structures, reference="label:chain")`の出力が完全に一致することを
テストで確認)。

## 出力セクションの絞り込み(FASTA・Substitutionsセクションの削除)

上記の対応の直後、ユーザーから「`== Sequences (FASTA, observed residues only) ==`セクションは
いりません。`== Substitutions relative to reference ==`セクションもいりません」との明確な指示を
受けた。レポートを「Pairwise identity」「Alignment(整列表示)」の2セクションのみに絞り込んだ。

- `build_report()`から`format_fasta()`/`format_mutation_report()`の呼び出しを削除。
- ただし`format_fasta()`・`format_mutation_report()`(および内部で使う`structcompare`/`seqalign`の
  置換・欠損検出ロジック)自体は削除せず残した(それぞれ単体テストが存在し、独立して再利用しうる
  ビルディングブロックであるため。BRAF V600E検出等ですでに実データ検証済みの機能でもあり、
  「セクションが不要」という指示を「機能自体の削除」まで拡大解釈しないよう留意した)。
- 一方、`--reference`のバリデーション(存在しないラベル/チェーンや不正な配列文字列でのエラー)は
  従来`format_mutation_report()`内で行っていたため、これを呼ばなくなると静かに握りつぶされて
  しまう(`format_alignment_block()`はreference行を追加できるかを静かに判定するだけで、
  無効な入力に対してエラーを出さない設計のため)。これを避けるため`_validate_reference()`を
  新設し、`build_report()`が`reference`指定時に必ず呼ぶことで、以前と同じエラーメッセージ
  (`chain not found: ...`/`--reference must be either ...`)を維持した。
- `--help`の説明文(`src/sequencealign/cli.py`)も、`--reference`の役割が「整列表示への基準行
  追加」に変わったことに合わせて全面的に書き直した。

## `--indir`解決ロジックの共通化

従来`alignview/cli.py`にのみ実装されていた`--indir`解決(拡張子省略時の自動補完、繰り返し指定、
"/"を含むトークンの絶対/相対パス扱い)を、2機能目(`sequencealign`)で必要になったタイミングで
`src/structio/resolve.py`(`resolve_structure_tokens()`)へ切り出し、`alignview/cli.py`を
書き換えて共用するようにした(既存テストは`CliRunner`経由でCLI全体の挙動を検証しているため
挙動に変更はない。専用ユニットテストを`tests/structio/test_resolve.py`に追加)。

## FASTA入力対応と`--reference`オプションの廃止

ユーザーからの要望: 「`pf sequence-align --indir data/mpro P0DTD1.fasta 6LU7_abc`のように、
FASTA形式も入力に対応してほしい。拡張子がない場合は優先順位`.cif`, `.pdb`, `.fasta`の順。
`P0DTD1`は3次元構造を持たずFASTAファイルのみのため、拡張子省略(`P0DTD1`のみ)でも機能してほしい」。

- `src/structio/resolve.py`の`resolve_structure_tokens()`に`extensions`引数を追加し、
  拡張子省略時に試す拡張子リストを呼び出し元ごとに変更できるようにした(既定は従来通り
  `.cif`/`.mmcif`/`.pdb`)。`sequencealign/cli.py`は`(".cif", ".mmcif", ".pdb", ".fasta")`を渡す。
  `alignview`(PyMOLでの3次元表示が前提)は`.fasta`を追加しておらず、従来通り構造ファイルの
  拡張子のみ解決する(共有インフラである`resolve.py`本体に`.fasta`をグローバルに追加すると、
  align-viewが構造を持たないFASTAを誤って解決してしまうため、呼び出し元ごとに切り替える設計にした)。
- `src/structio/fasta.py`を新設し、`parse_fasta(path) -> list[tuple[str, str]]`(Biopython`SeqIO`で
  ヘッダー・配列のタプル列を返す)をアトミックな技術要素として追加した。
- `LabeledStructure.atoms`を`Atomic | None`に変更し、`.fasta`から読み込んだ場合は`atoms=None`とした
  (3次元構造を持たないため)。`chains`はFASTAの各レコードをA, B, C...と順にチェーンIDを振った
  `ChainSequence`に変換して構築する(1残基目をresnum=1として連番。従来の`--reference`アミノ酸配列
  指定時と同じ前提)。
- `atoms`を要する処理(`format_identity_matrix()`の`matchChains`呼び出し、
  `format_mutation_report()`の`label:chain_id`基準)は、`atoms is None`の構造を対象から除外する
  (前者は静かにスキップ、後者は基準に指定された場合は明確な`ValueError`、比較対象に含まれる場合は
  「no atomic structure (loaded from FASTA)」という行を出す)。

この対応により、基準配列(UniProt正規配列等)を渡す専用の`--reference`オプションが不要になった
(整列表示への基準行追加という`--reference`の主用途を、`.fasta`を通常の入力ファイルとして渡す方式が
完全に代替するため)、とユーザーから指摘があり、`--reference`オプションを廃止した:

- `sequence_align_cmd`から`--reference`の`click.option`を削除。
- `build_report()`/`format_alignment_block()`から`reference`引数と、それに伴うreference行挿入ロジック・
  `_validate_reference()`を削除した(いずれも`--reference`専用のロジックで、他から呼ばれていなかった)。
- `format_mutation_report()`(`label:chain_id`/アミノ酸配列のどちらでも基準にできる置換検出)は
  `--reference`とは独立した別機能であり、そもそも`build_report()`(実際のCLI出力)からは呼ばれて
  いなかったため(上記「出力セクションの絞り込み」参照)、今回の廃止対象外とし変更していない
  (`atoms is None`のガードのみ追加)。

`--indir`解決の優先順位はユーザー指定通り`.cif` → `.pdb` → `.fasta`(既存の`.mmcif`は`.cif`と`.pdb`の
間に維持、実質的に「`.cif`系を最優先、次に`.pdb`、最後に配列のみの`.fasta`」という意図を保つ)。

## 実装ファイル

- `src/seqextract/chains.py` — 構造(Atomic)からの蛋白チェーン配列+残基番号の抽出(ProDy)
- `src/structcompare/compare.py` — 構造間のチェーン単位配列比較(ProDy `matchChains`ラッパー)
- `src/seqalign/pairwise.py` — 任意配列同士のペアワイズアラインメント(Biopython `PairwiseAligner`直接利用)
- `src/structio/resolve.py` — `--indir`解決ロジック(`alignview`と共用、対応拡張子は呼び出し元が指定)
- `src/structio/fasta.py` — FASTAファイルの読み込み(Biopython `SeqIO`、アトミックな技術要素として分離)
- `src/sequencealign/report.py` — レポート組み立て(sequencealign固有のドメインロジック)
- `src/sequencealign/cli.py` — `pf sequence-align` サブコマンド

## 依存パッケージの追加

`biopython`を`pharmoforge`環境に追加(`mamba install -n pharmoforge -c conda-forge biopython`)、
`pyproject.toml`の`dependencies`にも追記。ProDy自体がBiopythonをオプショナルに要求する
(`matchChains`等の内部で使用)ことと、`seqalign`での直接利用の両方の理由による。

## テスト

```bash
pytest tests/sequencealign tests/seqextract tests/structcompare tests/seqalign tests/structio tests/alignview
```

## 動作例(実データ)

```bash
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160
pf sequence-align --indir data/mpro P0DTD1.fasta 6LU7_abc
pf sequence-align --indir data/mpro P0DTD1 6LU7_abc
# => 上記2つは同じ結果になる(P0DTD1.cif/.pdbが存在しないため.fastaに解決される)。
#    P0DTD1:AがAlignmentセクションに他の構造と同じ行として加わる(reference専用の仕組みは廃止)。
```

(`--reference`オプション廃止前の動作例。BRAF(P15056)のAlphaFoldモデルと結晶構造4MNFを
`--reference P15056_AF:A`で比較すると既知の発がん性変異V600Eが正しく検出されることを確認した
記録は、上記「`--reference`オプションの仕様」「実データでの検証結果」節参照。現在のCLIには
`--reference`はないが、`format_mutation_report()`自体はライブラリ関数として残っており
`pytest tests/sequencealign -k format_mutation_report`で動作確認できる。)

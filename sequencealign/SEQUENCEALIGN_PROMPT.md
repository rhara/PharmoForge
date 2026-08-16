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

```
pf sequence-align <構造ファイル1> [<構造ファイル2> ...] [--reference <基準>] [--width <残基数>] [--output <出力ファイル>]
```

- 構造ファイルの指定は`align-view`と同じトークン解決(`--indir`、拡張子省略時の自動補完)を
  共有する(下記「`--indir`解決ロジックの共通化」参照)。
- `--reference`省略時は置換一覧セクションを出力しない。
- `--width`(既定`DEFAULT_ALIGN_WIDTH=100`)は整列表示セクションの折り返し残基数。
- `--output`/`-o`省略時は標準出力にレポートを出す。

### 出力レポートの構成(`src/sequencealign/report.py`)

レポート本文はすべて英語(下記「コマンド出力の英語化」参照)。

1. `== Sequences (FASTA, observed residues only) ==`: 各構造・各蛋白チェーンの配列
   (`seqextract.get_chain_sequences()`、CA原子ベース)。
2. `== Pairwise identity ==`: 全構造の組み合わせについて`structcompare.match_chains()`の結果を一覧化。
3. `== Alignment (by residue number) ==`: 残基番号ベースの整列表示(上記「整列表示」参照)。
4. `== Substitutions relative to reference ==`(`--reference`指定時のみ): 見出し行
   (`reference: ...`/`reference sequence: ...`)に続けて基準配列自体をFASTA形式(60残基/行)で
   出力した上で、上記の判定に応じて`structcompare.find_substitutions()`または
   `seqalign.align_to_reference()`を使用した置換一覧を出力する。置換一覧の各行に続けて、
   基準に対する欠損領域(「gaps:」行、上記「欠損領域(ギャップ)の可視化」参照)を出力する。
   基準配列自体の出力は、ユーザーからの「referenceのシーケンスも出力してほしい」という
   要望を受けて追加した(従来は`ラベル:チェーンID`基準の場合はFASTAセクションを見れば
   分かったが、アミノ酸配列を直接指定した場合はその配列自体がレポートのどこにも
   出力されていなかった)。

## `--indir`解決ロジックの共通化

従来`alignview/cli.py`にのみ実装されていた`--indir`解決(拡張子省略時の自動補完、繰り返し指定、
"/"を含むトークンの絶対/相対パス扱い)を、2機能目(`sequencealign`)で必要になったタイミングで
`src/structio/resolve.py`(`resolve_structure_tokens()`)へ切り出し、`alignview/cli.py`を
書き換えて共用するようにした(既存テストは`CliRunner`経由でCLI全体の挙動を検証しているため
挙動に変更はない。専用ユニットテストを`tests/structio/test_resolve.py`に追加)。

## 実装ファイル

- `src/seqextract/chains.py` — 構造(Atomic)からの蛋白チェーン配列+残基番号の抽出(ProDy)
- `src/structcompare/compare.py` — 構造間のチェーン単位配列比較(ProDy `matchChains`ラッパー)
- `src/seqalign/pairwise.py` — 任意配列同士のペアワイズアラインメント(Biopython `PairwiseAligner`直接利用)
- `src/structio/resolve.py` — `--indir`解決ロジック(`alignview`と共用)
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
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a --reference P24941_AF:A
# => 1AQ1_ab: no substitutions (...)
#      gaps: reference only (missing in target): 36-43, 149-161
pf sequence-align --indir data/braf P15056_AF 4MNF_ac --reference P15056_AF:A
# => 4MNF_ac: 1 substitution(s) (seqid=99.6%, overlap=33.6%): V600E
#      gaps: reference only (missing in target): 1-448, 601-615, 721-766
pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab 1HCL_a
pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160
```

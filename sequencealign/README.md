# sequencealign

複数のPDB/CIF構造から蛋白チェーンの配列を抽出し、FASTA出力・pairwise配列同一性(%identity)・
残基番号ベースの整列表示・基準配列に対する残基置換(変異)一覧を出力する機能。ダウンロードした
複数の構造(結晶構造・AlphaFold予測構造)間で、配列が本当に同一か、どこに変異・構築上の違いが
あるかを素早く確認するためのもの。

## 使い方

```bash
pf sequence-align <構造ファイル1> <構造ファイル2> [...] [--reference <基準>] [--width <残基数>] [--output <出力ファイル>]
```

```bash
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a --reference P24941_AF:A
pf sequence-align data/braf/P15056_AF.cif data/braf/3OG7_ac.cif --reference P15056_AF:A -o report.txt
pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab 1HCL_a
pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160
```

コマンドの出力(レポート本文・`--help`)はすべて英語。README・PROMPT等のドキュメントは日本語。

- 構造ファイルはファイル名(拡張子抜き)をラベルとして扱う(FASTAヘッダ・identity表・置換一覧で
  `<ラベル>:<チェーンID>`として参照する)。
- `--indir DIR`: `align-view`と同じ引数体系。繰り返し指定可能で、以降のファイル名(拡張子省略可、
  `.cif`優先、次に`.pdb`)を`DIR`配下から解決する。`/`を含む指定(または絶対パス)は`--indir`に
  よらずカレントディレクトリ相対 or 絶対パスとして扱う。
- `--width`(既定`100`): 整列表示セクションの折り返し残基数。
- `--output` / `-o`: レポートの保存先(省略時は標準出力)。
- `--reference`: 指定すると残基置換一覧を出力する(省略時は出力しない)。2通りの指定方法:
  - `ラベル:チェーンID`(例: `P24941_AF:A`): 読み込んだ構造の1チェーンを基準にする。
    残基番号ベースの対応付けのみを用いる(同一蛋白の構造間ではPDBの残基番号が揃っている前提。
    [`pf align-view --method number`](../alignview/README.md)と同じ前提)。番号体系が異なる
    構造間では「対応が取れませんでした」と表示される。
  - アミノ酸配列(1文字表記、コロンを含まない文字列): 構造を伴わない任意配列(UniProt正規配列や
    ユーザー指定の基準配列)を直接指定できる。この場合は配列アラインメント([`seqalign`](../API.md#srcseqalign))
    を用いるため、残基番号が揃っていない構造間でも比較できる。

### 出力

4つのセクションからなるテキストレポートを出力する。

1. **配列(FASTA)**: 各構造・各蛋白チェーンの配列(観測されたCA原子のみ、電子密度が見えず欠損した
   残基は含まれない。UniProtの完全配列とは異なりうる)。ヘッダは`><ラベル>:<チェーンID>
   length=<配列長> range=<開始残基番号>-<終了残基番号>`。
2. **Pairwise identity**: 全構造の組み合わせについて、チェーン単位の%identity/%overlapを一覧化
   ([`structcompare`](../API.md#srcstructcompare)、`matchChains`)。
3. **整列表示(残基番号ベース)**: 全構造・全蛋白チェーンの配列を、残基番号を共通の軸として
   縦に並べて表示する(`--width`残基/行、既定100で折り返し)。各ブロックの直上に10残基ごとの
   位置番号+`|`の目盛り(ルーラー)を表示する(右揃えでは収まらない目盛りは、`|`は本来の列の
   まま数字だけをブロック左端に寄せて表示し、値の欠落・誤読を防ぐ)。配列アラインメントは行わず、
   構造間でPDBの残基番号が揃っている前提で並べる(`pf align-view --method number`と同じ前提)。
   観測されていない残基は`-`で埋める。異なる蛋白の構造を混在させると無意味な結果になるため、
   通常は同一蛋白の複数構造を
   対象とする。
4. **基準配列に対する置換**(`--reference`指定時のみ): 基準に対する各構造の残基置換
   (例: `V600E`)一覧。`ラベル:チェーンID`基準の場合は`資源:残基番号:置換後`、配列基準の場合は
   基準配列内の位置と構造側の実際の残基番号を併記する(例: `V600E(構造残基番号=600)`)。
   置換とは別に、基準に対して欠損している(または基準にない)領域を`欠損:`行で示す
   (PDB構造は電子密度が見えない領域(N/C末端の未解析部分・ループ等)が欠落することが多いため)。

## 実装方針

- `src/seqextract`(配列+残基番号の抽出、ProDy)・`src/structcompare`(構造間の配列比較、ProDyの
  `matchChains`ラッパー)・`src/seqalign`(任意配列同士の比較、Biopython`PairwiseAligner`)という
  3つのアトミックなパッケージを組み合わせて`src/sequencealign/report.py`でレポートを組み立てる
  (詳細は[API.md](../API.md)参照)。
- `--indir`解決ロジックは[`structio.resolve`](../API.md#srcstructio)を`align-view`と共用する。

### `pwalign`に関する既知の注意点(`structcompare`)

ProDyの`matchChains`は残基番号ベースの直接対応付け(`pwalign=False`)を基本とし、対応が取れない
場合にBiopythonによる配列アラインメントへフォールバックする(`pwalign=True`)。現行のprodyバージョン
(2.6.1)では、`pwalign=True`時に返る`AtomMap`の残基番号が実際の対応関係を反映しない(重ね合わせ後
RMSDが数十Åに達する)ことを実データで確認しているため、**残基単位の対応(置換検出)には
`pwalign=False`のみを用いる**。`pwalign=True`は%identity/%overlapの把握(pairwise identity表)
にのみ用いる(この値自体は`pwalign`の有無によらず一致することを確認済み)。また`pwalign=True`は
入力の組み合わせによって内部で`StopIteration`を送出して失敗することがあるため、その場合は
警告ログを出しつつ「対応なし」として扱う。

## テスト

```bash
pytest tests/sequencealign tests/seqextract tests/structcompare tests/seqalign tests/structio
```

## 動作例(実データ、動作確認済み)

BRAF(P15056)のAlphaFoldモデルと結晶構造4MNF/3OG7を比較すると、既知の発がん性変異V600Eが
正しく検出され、かつ結晶構造でN/C末端(全長のうちキナーゼドメイン以外)や不可視のループが
欠損していることも合わせてわかる:

```bash
pf sequence-align --indir data/braf P15056_AF 4MNF_ac 3OG7_ac --reference P15056_AF:A
# ...
# == Substitutions relative to reference ==
# reference: P15056_AF:A
#   4MNF_ac: 1 substitution(s) (seqid=99.6%, overlap=33.6%): V600E
#     gaps: reference only (missing in target): 1-448, 601-615, 721-766
#   3OG7_ac: 14 substitution(s) (seqid=94.3%, overlap=32.2%): K522A, I543A, ...
#     gaps: reference only (missing in target): 1-448, 545-547, 597-614, 627-630, 721-766
```

CDK2(P24941)のAlphaFoldモデルと結晶構造2件の整列表示(残基番号ベース、ルーラー付き)の例:

```bash
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
# == Alignment (by residue number) ==
#                      10        20        30        40        50        60        70        80        90       100
#                       |         |         |         |         |         |         |         |         |         |
# P24941_AF:A  MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
# 1AQ1_ab:A    MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKI--------VPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
# 1HCL_a:A     MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIR----TEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
```

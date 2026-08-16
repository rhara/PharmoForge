# sequencealign

複数のPDB/CIF構造から蛋白チェーンの配列を抽出し、pairwise配列同一性(%identity)・
残基番号ベースの整列表示を出力する機能。ダウンロードした複数の構造(結晶構造・AlphaFold
予測構造)間で、配列が本当に同一か、どこに構築上の違いがあるかを素早く確認するためのもの。

## 使い方

```bash
pf sequence-align <構造ファイル1> <構造ファイル2> [...] [--reference <基準>] [--width <残基数>] [--output <出力ファイル>]
```

```bash
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab 1HCL_a
pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160 -o report.txt
```

コマンドの出力(レポート本文・`--help`)はすべて英語。README・PROMPT等のドキュメントは日本語。

- 構造ファイルはファイル名(拡張子抜き)をラベルとして扱う(identity表・整列表示で
  `<ラベル>:<チェーンID>`として参照する)。
- `--indir DIR`: `align-view`と同じ引数体系。繰り返し指定可能で、以降のファイル名(拡張子省略可、
  `.cif`優先、次に`.pdb`)を`DIR`配下から解決する。`/`を含む指定(または絶対パス)は`--indir`に
  よらずカレントディレクトリ相対 or 絶対パスとして扱う。
- `--width`(既定`100`): 整列表示セクションの折り返し残基数。
- `--output` / `-o`: レポートの保存先(省略時は標準出力)。
- `--reference`: 整列表示に基準配列の行を加える。2通りの指定方法:
  - `ラベル:チェーンID`(例: `P24941_AF:A`): 読み込んだ構造の1チェーン。そのチェーンは
    すでに構造側の行として整列表示に含まれているため、存在確認(存在しなければエラー)のみ
    行い、行は追加しない。
  - アミノ酸配列(1文字表記、コロンを含まない文字列): 構造を伴わない任意配列(UniProt正規配列や
    ユーザー指定の基準配列)を`reference`という行として整列表示に追加する。配列の1文字目を
    残基番号1として扱う(基準配列が構造と同じ番号体系で、通常は残基1から始まる前提。
    [`pf align-view --method number`](../alignview/README.md)と同じ前提)。

### 出力

2つのセクションからなるテキストレポートを出力する。

1. **Pairwise identity**: 全構造の組み合わせについて、チェーン単位の%identity/%overlapを一覧化
   ([`structcompare`](../API.md#srcstructcompare)、`matchChains`)。
2. **整列表示(残基番号ベース)**: 全構造・全蛋白チェーンの配列(観測されたCA原子のみ、電子密度が
   見えず欠損した残基は`-`で埋める)を、残基番号を共通の軸として縦に並べて表示する(`--width`
   残基/行、既定100で折り返し)。各ブロックの直上に10残基ごとの位置番号+`|`の目盛り(ルーラー)を
   表示する(右揃えでは収まらない目盛りは、`|`は本来の列のまま数字だけをブロック左端に寄せて
   表示し、値の欠落・誤読を防ぐ)。配列アラインメントは行わず、構造間でPDBの残基番号が揃っている
   前提で並べる。異なる蛋白の構造を混在させると無意味な結果になるため、通常は同一蛋白の複数構造を
   対象とする。`--reference`にアミノ酸配列を直接指定した場合は、`reference`という行を追加する。

## 実装方針

- `src/seqextract`(配列+残基番号の抽出、ProDy)・`src/structcompare`(構造間の配列比較、ProDyの
  `matchChains`ラッパー)・`src/seqalign`(任意配列同士の比較、Biopython`PairwiseAligner`)という
  3つのアトミックなパッケージを組み合わせて`src/sequencealign/report.py`でレポートを組み立てる
  (詳細は[API.md](../API.md)参照)。`structcompare.match_chains()`はPairwise identityセクションで
  引き続き使う。一方`structcompare.find_substitutions()`/`seqalign.align_to_reference()`(置換・
  欠損検出)は現在レポートからは使われていないが(下記「出力セクションの絞り込み」参照)、
  関数自体は残しておりテストもある。
- `--indir`解決ロジックは[`structio.resolve`](../API.md#srcstructio)を`align-view`と共用する。

### 出力セクションの絞り込み

当初はFASTA配列一覧・基準配列に対する残基置換(変異)/欠損領域一覧のセクションもあったが、
ユーザーから「`== Sequences ==`と`== Substitutions relative to reference ==`のセクションは
いらない」との要望を受けて出力から外した(詳細は[SEQUENCEALIGN_PROMPT.md](SEQUENCEALIGN_PROMPT.md)
参照)。`--reference`は現在、整列表示に基準配列の行を加える目的にのみ使う。

## テスト

```bash
pytest tests/sequencealign tests/seqextract tests/structcompare tests/seqalign tests/structio
```

## 動作例(実データ、動作確認済み)

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

`--reference`にアミノ酸配列を直接指定した場合は、`reference`行としてこの整列表示にも加わる:

```bash
pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab
# == Alignment (by residue number) ==
#                    10        20        30        40        50        60        70        80        90       100
#                     |         |         |         |         |         |         |         |         |         |
# reference  MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
# 1AQ1_ab:A  MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKI--------VPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
```

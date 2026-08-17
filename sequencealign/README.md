# sequencealign

複数のPDB/CIF構造(および/またはFASTA配列)から蛋白チェーンの配列を抽出し、pairwise配列同一性
(%identity/%coverage)・配列アラインメントに基づく整列表示を出力する機能。ダウンロードした複数の
構造(結晶構造・AlphaFold予測構造)間で、配列が本当に同一か、どこに構築上の違いがあるかを
素早く確認するためのもの。

## 使い方

```bash
pf sequence-align <構造/FASTAファイル1> <構造/FASTAファイル2> [...] [--method align|number] [--identity-format combined|separate] [--width <残基数>] [--output <出力ファイル>]
```

```bash
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160 -o report.txt
pf sequence-align --indir data/mpro P0DTD1.fasta 6LU7_abc
pf sequence-align --indir data/mpro P0DTD1 6LU7_abc
pf sequence-align --indir data/mpro P0DTD1 6LU7_abc --identity-format separate
```

コマンドの出力(レポート本文・`--help`)はすべて英語。README・PROMPT等のドキュメントは日本語。

- 構造/FASTAファイルはファイル名(拡張子抜き)をラベルとして扱う(identity表・整列表示で
  `<ラベル>:<チェーンID>`として参照する)。
- `--indir DIR`: `align-view`と同じ引数体系。繰り返し指定可能で、以降のファイル名(拡張子省略可、
  `.cif`→`.mmcif`→`.pdb`→`.fasta`の順で解決)を`DIR`配下から解決する。`/`を含む指定(または絶対パス)は
  `--indir`によらずカレントディレクトリ相対 or 絶対パスとして扱う。
- `--method`(`align`/`number`、既定`align`): 整列表示セクションの並べ方。下記参照。
- `--identity-format`(`combined`/`separate`、既定`combined`): Pairwise identityセクションの
  表示形式。下記参照。
- `--width`(既定`100`): 整列表示セクションの折り返し残基数。
- `--output` / `-o`: レポートの保存先(省略時は標準出力)。
- FASTA入力(`.fasta`): 3次元構造を持たないが、Pairwise identity・整列表示のどちらのセクションにも
  他の構造と同じ行として含まれる(1残基目をresnum=1として連番を振る)。基準配列(UniProt正規配列等)を
  加えたい場合は、その配列を`.fasta`として他の構造ファイルと同様に入力に含めればよい(専用の
  `--reference`オプションは廃止した。詳細は下記「出力セクションの絞り込み」参照)。1つのFASTAに
  複数レコードがある場合は、A, B, C...と順にチェーンIDを振る。

### 出力

`--identity-format`により1〜2個のPairwise identity/coverage系セクション+1個の整列表示セクションを
出力する。

1. **Pairwise identity(/coverage)**: 全チェーンの組み合わせ(構造・FASTA・同一構造内の複数チェーン
   同士を問わず全て)について、ペアワイズグローバルアラインメント([`seqalign.align_to_reference`]
   (../API.md#srcseqalign)、Biopython `PairwiseAligner`)による%identity・%coverageをN×Nの
   グリッド表として一覧化する(対角は`-`)。
   - `--identity-format combined`(既定): 1つの表(`== Pairwise identity/coverage ==`)に
     `identity/coverage`(%記号なし)をセルとして表示する。
   - `--identity-format separate`: `== Pairwise identity ==`(identityのみ)と`== Coverage ==`
     (coverageのみ)の2つの表に分ける。
   - coverageは**非対称**な値である点に注意: セル`(行, 列)`は「行のチェーンを基準配列としたとき、
     アラインメントで列のチェーンと対応した割合」(基準配列=行の長さが分母)。例えば全長配列(行)を
     ドメインのみの構造(列)と比較すると低coverageに、逆にドメインのみの構造(行)を全長配列(列)と
     比較すると高coverageになる。identityはほぼ対称(どちらを基準にしてもアラインメントされた
     位置の一致率はほぼ同じ)なので非対称の扱いはしていない。
2. **整列表示**: 全構造・全蛋白チェーンの配列(構造はCA原子のみ観測、電子密度が見えず欠損した残基は
   `-`で埋める。FASTAはそのままの配列)を、共通の軸を基準に縦に並べて表示する(`--width`残基/行、
   既定100で折り返し)。各ブロックの直上に10残基ごとの位置番号+`|`の目盛り(ルーラー)を表示する
   (右揃えでは収まらない目盛りは、`|`は本来の列のまま数字だけをブロック左端に寄せて表示し、値の
   欠落・誤読を防ぐ)。異なる蛋白の配列を混在させると無意味な結果になるため、通常は同一蛋白の
   複数構造・配列を対象とする。軸の決め方は`--method`で選ぶ:
   - `align`(既定、`== Alignment (sequence-aligned) ==`): 先頭の入力の先頭チェーンを基準に、
     他の全チェーンをペアワイズグローバルアラインメント(`seqalign.align_to_reference`)した結果で
     位置を揃える。構造間でPDBの残基番号が揃っていなくても(例: 全長配列 vs ドメインのみの構造)
     正しい位置に並ぶ。複数配列同時アラインメント(MSA)ではなく、先頭チェーンに対する個別の
     ペアワイズアラインメントである点に注意。
   - `number`(`== Alignment (by residue number) ==`): 配列アラインメントを行わず、残基番号が
     一致する列に同じアミノ酸が並ぶ前提で並べる(構造間でPDBの残基番号が既に揃っている、例えば
     同じUniProt番号体系であることが分かっている場合のみ有効)。

## 実装方針

- `src/seqextract`(配列+残基番号の抽出、ProDy)・`src/seqalign`(任意配列同士のペアワイズ
  グローバルアラインメント、Biopython`PairwiseAligner`)という2つのアトミックなパッケージを
  組み合わせて`src/sequencealign/report.py`でレポートを組み立てる(詳細は[API.md](../API.md)参照)。
  Pairwise identity/coverage・整列表示(`--method align`)のいずれも`seqalign.align_to_reference()`
  を使う(`structcompare.match_chains()`は現在使っていない)。`structcompare.find_substitutions()`/
  `structcompare.match_chains()`(構造=ProDy Atomic同士の比較)は現在レポートの既定出力からは
  使われていないが(下記「出力セクションの絞り込み」参照)、関数自体は残しておりテストもある。
- `--indir`解決ロジックは[`structio.resolve`](../API.md#srcstructio)を`align-view`と共用する
  (対応拡張子は`sequence-align`側で`.fasta`を追加した独自リストを渡す)。
- FASTAの読み込みは[`structio.parse_fasta`](../API.md#srcstructio)(Biopython`SeqIO`)を使う。

### 出力セクションの絞り込み・`--reference`の廃止

当初はFASTA配列一覧・基準配列に対する残基置換(変異)/欠損領域一覧のセクションもあったが、
ユーザーから「`== Sequences ==`と`== Substitutions relative to reference ==`のセクションは
いらない」との要望を受けて出力から外した(詳細は[SEQUENCEALIGN_PROMPT.md](SEQUENCEALIGN_PROMPT.md)
参照)。その後、入力ファイルとして`.fasta`を直接渡せるようになったことで、基準配列を渡す専用の
`--reference`オプション(整列表示に`reference`行を加える機能)は不要になったため廃止した。
基準配列を加えたい場合は、他の構造ファイルと同様にその配列を`.fasta`として入力トークンに含めれば、
`<ラベル>:A`という通常の行として同じ整列表示に加わる。

### 残基番号ベース整列の限界と`--method align`・Pairwise identityのグリッド化

当初、整列表示は残基番号ベース(現在の`--method number`)のみだった。FASTA入力対応後、構造間で
残基番号体系が異なる組み合わせ(例: 全長のUniProt配列 vs 一部ドメインのみの結晶構造)では、
残基番号ベースの整列が意味をなさない(無関係な位置同士が同じ列に並んでしまう)ことが実データで
判明した。これを受けて`seqalign.align_to_reference()`による配列アラインメントベースの整列
(`--method align`)を追加し、既定にした。Pairwise identityセクションも、当初は`structcompare`
(構造=Atomic同士の比較)ベースで、atomsを持たないFASTA入力を除外していたが、`--method align`と
同じ`seqalign.align_to_reference()`を使う方式に統一し、FASTA入力も含めた全チェーンのN×Nグリッド
(基準構造も含む)として一覧化するよう変更した。identity/coverageの2軸をどう表示するか
(1表に統合`combined` or 2表に分離`separate`)はユーザーの好みが定まらなかったため、両方式を
`--identity-format`オプションで選べるようにして残している。詳細は
[SEQUENCEALIGN_PROMPT.md](SEQUENCEALIGN_PROMPT.md)参照。

## テスト

```bash
pytest tests/sequencealign tests/seqextract tests/structcompare tests/seqalign tests/structio
```

## 動作例(実データ、動作確認済み)

SARS-CoV-2のreplicase polyprotein 1ab全長(P0DTD1、UniProt由来FASTA)とMproドメインのみの
結晶構造(6LU7)を組み合わせた例(`--method align`(既定)でないと正しく整列しない):

```bash
pf sequence-align --indir data/mpro P0DTD1 6LU7
# == Pairwise identity/coverage ==
#               P0DTD1:A      6LU7:A      6LU7:C
# P0DTD1:A             -   100.0/4.3   100.0/0.0
# 6LU7:A     100.0/100.0           -   100.0/1.0
# 6LU7:C     100.0/100.0 100.0/100.0           -
#
# == Alignment (sequence-aligned) ==
# (中略、6LU7:Aの配列が現れ始めるブロック)
#                3310      3320      3330      3340      3350      3360      3370      3380      3390      3400
#                   |         |         |         |         |         |         |         |         |         |
# P0DTD1:A  CPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIK
# 6LU7:A    CPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIK
```

`6LU7`(Mproドメインのみの結晶構造、ローカル採番1-306)が、全長`P0DTD1`(7096残基)のうち
実際にMproが位置する残基3264付近から正しく整列している(既知の生物学的事実と一致。
残基番号ベース(`--method number`)ではこの組み合わせは正しく整列できない)。

CDK2(P24941)のAlphaFoldモデルと結晶構造2件の整列表示(残基番号体系が揃っている場合、
`--method number`でも動く例):

```bash
pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a --method number
# == Alignment (by residue number) ==
#                      10        20        30        40        50        60        70        80        90       100
#                       |         |         |         |         |         |         |         |         |         |
# P24941_AF:A  MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
# 1AQ1_ab:A    MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKI--------VPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
# 1HCL_a:A     MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIR----TEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIP
```

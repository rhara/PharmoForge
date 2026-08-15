# scaffoldanalyzer

活性データ(化合物SMILES + 活性値)を、Bemis-Murckoスキャフォールド単位でクラスタリングし、
高活性群・低活性群それぞれにどのスキャフォールドが偏って出現するかを比較する機能。
`pf fetch activities=...` の出力をそのまま入力にできる。

## 使い方

```bash
pf analyze-scaffolds <入力TSV/CSV> --output-dir <出力ディレクトリ>
```

```bash
pf analyze-scaffolds data/cdk4_human_activities.tsv --output-dir data/cdk4_scaffold_analysis
```

### 処理内容

1. 各化合物のSMILES(既定で`smiles`列、`--smiles-col`で変更可)からBemis-Murckoスキャフォールドを求める。
   パースできない行は除外する(除外数はログに表示)。
2. 活性値列(既定`_median`、`--activity-col`で変更可)の分位点で各化合物を`high`/`mid`/`low`に分類する
   (既定: 上位25%を`high`、下位25%を`low`、`--high-quantile`/`--low-quantile`で変更可)。
3. スキャフォールドごとにhigh/mid/low件数と活性値のmean/medianを集計し、
   `enrichment = (スキャフォールドのhigh群内での出現割合) - (low群内での出現割合)` を計算する。
   出現総数が`--min-count`(既定2)未満のスキャフォールドは、解釈の信頼性が低いため集計から除外する。
4. `enrichment`降順で全スキャフォールドをTSVに出力し、上位・下位`--top-n`件(既定20)を
   構造式グリッド画像として出力する。
5. 同じ上位・下位`--top-n`件のスキャフォールドについて、それぞれに属する個々の化合物
   (置換基込み)を構造式付きのHTMLテーブルとして出力する。スキャフォールドごとに見出し行
   (構造式・n_total等の統計)を挟み、グループ内は活性値降順で化合物を並べる。

### 出力

`--output-dir`配下に以下を出力する。

| ファイル | 内容 |
| --- | --- |
| `scaffold_summary.tsv` | 全スキャフォールドの集計結果(`scaffold, n_total, n_high, n_mid, n_low, frac_high, frac_low, enrichment, mean_activity, median_activity`)、`enrichment`降順 |
| `scaffold_grid_high.png` | `enrichment`上位(高活性群に偏って出現)のスキャフォールドの構造式グリッド |
| `scaffold_grid_low.png` | `enrichment`下位(低活性群に偏って出現)のスキャフォールドの構造式グリッド |
| `scaffold_compounds_high.html` | `enrichment`上位スキャフォールドに属する個々の化合物(置換基込み)の構造式・SMILES・活性値・binの一覧表 |
| `scaffold_compounds_low.html` | `enrichment`下位スキャフォールドに属する個々の化合物の一覧表 |

## 実装方針

- スキャフォールド計算は共通パッケージ[`src/molscaffold`](../src/molscaffold)で行う(RDKitの`MurckoScaffold`のみ使用、側鎖除去・環の飽和度等はそのまま)。アトミックな技術要素として独立させている。
- 活性値の分位点ビニングは共通パッケージ[`src/actbin`](../src/actbin)で行う。同じくアトミックな技術要素として独立させている。
- スキャフォールド単位での集計(high/mid/low件数、enrichment計算)は本機能固有のロジック(`src/scaffoldanalyzer/clustering.py`)。
- グリッド画像・化合物構造画像は`rdkit.Chem.Draw`(`MolsToGridImage`/`MolToImage`、内部で`pillow`を使用)。化合物テーブルの構造式はPNGをbase64でHTMLに埋め込み、外部ファイル依存なしの単一HTMLとして出力する。

共通パッケージ側の関数一覧は[API.md](../API.md)を参照。

## テスト

```bash
pytest tests/scaffoldanalyzer tests/molscaffold tests/actbin
```

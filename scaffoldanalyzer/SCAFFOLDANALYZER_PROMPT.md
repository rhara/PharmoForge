# scaffoldanalyzer 実装記録

このドキュメントは `scaffoldanalyzer` 機能を再現するための仕様記録。

## 目的

活性データ(化合物SMILES + 活性値)を、Bemis-Murckoスキャフォールド単位でクラスタリングし、
高活性群・低活性群それぞれにどのスキャフォールドが偏って出現するかを比較する。
`pf fetch activities=...`(fetcher機能)の出力をそのまま入力にできる、探索的なSAR解析機能。

## CLI仕様

```
pf analyze-scaffolds <入力TSV/CSV> --output-dir <出力ディレクトリ>
```

```
pf analyze-scaffolds data/cdk4_human_activities.tsv --output-dir data/cdk4_scaffold_analysis
```

オプション:

| オプション | 既定値 | 説明 |
| --- | --- | --- |
| `--smiles-col` | `smiles` | SMILES列名 |
| `--activity-col` | `_median` | 活性値列名 |
| `--output-dir` / `-o` | (必須) | 出力ディレクトリ |
| `--high-quantile` | `0.75` | 高活性とみなす分位点(以上) |
| `--low-quantile` | `0.25` | 低活性とみなす分位点(以下) |
| `--min-count` | `2` | 集計対象とするスキャフォールドの最小出現数(全bin合計) |
| `--top-n` | `20` | グリッド画像に描画する上位・下位件数 |

## 処理内容

1. 入力(TSV/CSV、区切り文字は`pandas.read_csv(sep=None, engine="python")`で自動判定)の各化合物SMILESから
   Bemis-Murckoスキャフォールドを計算する(`molscaffold.compute_scaffold()`)。パース失敗行は除外しログに出す
   (`scaffoldanalyzer.clustering.add_scaffolds()`)。
2. 活性値列の`--high-quantile`/`--low-quantile`分位点で各化合物を`high`/`mid`/`low`に分類する
   (`actbin.assign_activity_bins()`)。
3. スキャフォールドごとにhigh/mid/low件数と活性値のmean/medianを集計し、
   `enrichment = (high群内での出現割合) - (low群内での出現割合)` を計算する
   (`scaffoldanalyzer.clustering.summarize_scaffolds()`)。出現総数が`--min-count`未満のスキャフォールドは除外。
4. `enrichment`降順で全件をTSVに出力(`scaffoldanalyzer.report.write_summary_tsv()`)、
   上位・下位`--top-n`件を構造式グリッド画像として出力する(`scaffoldanalyzer.report.render_scaffold_grid()`、
   `RDKit.Draw.MolsToGridImage`)。
5. 同じ上位・下位`--top-n`件のスキャフォールドについて、それに属する個々の化合物(置換基込み)を
   構造式付きHTMLテーブルとして出力する(`scaffoldanalyzer.report.render_compound_table()`)。
   スキャフォールドごとに見出し行(構造式・n_total等)を挟み、グループ内は活性値降順で化合物を並べる。
   構造式は`RDKit.Draw.MolToImage`でPNG化しbase64埋め込み(外部ファイル依存なしの単一HTML)。

## 出力

`--output-dir`配下:

| ファイル | 内容 |
| --- | --- |
| `scaffold_summary.tsv` | 全スキャフォールドの集計結果(`scaffold, n_total, n_high, n_mid, n_low, frac_high, frac_low, enrichment, mean_activity, median_activity`)、`enrichment`降順 |
| `scaffold_grid_high.png` | `enrichment`上位のスキャフォールドの構造式グリッド(n/high/low件数、enrichment、median活性を凡例表示) |
| `scaffold_grid_low.png` | `enrichment`下位のスキャフォールドの構造式グリッド |
| `scaffold_compounds_high.html` | `enrichment`上位スキャフォールドに属する個々の化合物一覧(構造式・SMILES・活性値・bin) |
| `scaffold_compounds_low.html` | `enrichment`下位スキャフォールドに属する個々の化合物一覧 |

## 実装ファイル

- `src/scaffoldanalyzer/clustering.py` — スキャフォールド付与・集計ロジック(本機能固有)
- `src/scaffoldanalyzer/report.py` — TSV/グリッド画像出力(本機能固有)
- `src/scaffoldanalyzer/cli.py` — `pf analyze-scaffolds` サブコマンド
- `src/molscaffold/scaffold.py` — Bemis-Murckoスキャフォールド計算(アトミックな技術要素として分離、`compute_scaffold()`)
- `src/actbin/binning.py` — 活性値の分位点ビニング(アトミックな技術要素として分離、`assign_activity_bins()`)

## テスト

```bash
pytest tests/scaffoldanalyzer tests/molscaffold tests/actbin
```

## 動作例(サンプルデータ)

CDK4(ヒト)の活性データ([fetcher](../fetcher/FETCHER_PROMPT.md)で取得したもの)を題材にした解析例:

```bash
pf analyze-scaffolds data/cdk4_human_activities.tsv --output-dir data/cdk4_scaffold_analysis
```

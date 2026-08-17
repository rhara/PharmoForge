# pocketfinder

[fpocket](https://github.com/Discngine/fpocket)により蛋白構造からリガンド結合ポケット候補を検出し、
各ポケットに面する残基を一覧化する機能。

## 使い方

```bash
pf find-pocket <構造ファイル> --outdir <出力ディレクトリ> [--top N]
```

```bash
pf find-pocket data/cdk2/P24941_AF.cif --outdir data/cdk2/P24941_AF_pockets
pf find-pocket data/cdk2/1HCL.cif --outdir data/cdk2/1HCL_pockets --top 3
```

入力構造は[`pf fetch structure=...`](../fetcher/README.md)(RCSB PDB/AlphaFold DB)や
[`pf protein-extract`](../proteinextract/README.md)の出力をそのまま使える(水分子・不要チェーンの
除去は本機能の対象外。事前に`pf protein-extract --remove-water`等で済ませておく)。

### 出力

`--outdir`配下に以下が生成される。

- `pockets.json` — 検出されたポケットの一覧(スコア降順)。各ポケットの主要記述子
  (`score`/`druggability_score`/`n_alpha_spheres`/`volume`/`total_sasa`/`polar_sasa`/`apolar_sasa`/
  `hydrophobicity_score`)と、ポケットに面する残基一覧(`residues`: `chain_id`/`resnum`/`resname`、
  重複除去・chain/resnum昇順)を含む。
- `<構造ファイル名>` — 入力構造のコピー(fpocketが入力と同じディレクトリに出力を書くための作業用)。
- `<構造名>_out/` — fpocketの生出力一式(PyMOL/VMDの可視化スクリプト、pocketごとのalphaスフィア座標
  (`.pqr`)等)。ポケットを3次元的に確認したい場合はここのスクリプトをそのまま利用できる。

`--top N`(省略時は全件): スコア上位N件のみ`pockets.json`に出力する。

## ポケットの可視化(`pf view-pocket`)

```bash
pf view-pocket <構造ファイル> --pockets <pockets.jsonのパス> [--top N] [--pymol-env <env名>]
```

```bash
pf view-pocket data/cdk2/P24941_AF.cif --pockets data/cdk2/P24941_AF_pockets/pockets.json
pf view-pocket data/cdk2/P24941_AF.cif --pockets data/cdk2/P24941_AF_pockets/pockets.json --top 3
```

`pf find-pocket`が出力した`pockets.json`を読み込み、対象の構造をPyMOLで開いてポケットごとに
周辺残基を配色(スティック表示+半透明サーフェス)して強調表示する。ポケットごとに異なる色が
自動で割り当てられ、強調表示した残基全体にズームする。`--top N`(省略時は全件): スコア上位N件の
ポケットのみ強調表示する。`--pymol-env`(既定`pymol`)は下記のPyMOL実行環境と共通。

## ポケット一覧の表形式出力(`pf show-pocket`)

```bash
pf show-pocket <pockets.jsonのパス> [--output <出力TSVファイル>]
```

```bash
pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json
pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json --output data/cdk2/P24941_AF_pockets/pockets.tsv
pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json | column -t
```

`pf find-pocket`が出力した`pockets.json`を、1行1ポケットのTSVに整形して出力する。列は
`pocket_id`/`score`/`druggability_score`/`n_alpha_spheres`/`volume`/`total_sasa`/`polar_sasa`/
`apolar_sasa`/`hydrophobicity_score`/`n_residues`/`residues`(周辺残基を`<chain>:<resnum>`形式で
カンマ結合した1セル)。`--output`省略時は標準出力にTSVを出力する(`column -t`等へのパイプ利用を想定)。

## PyMOL/fpocketの実行環境

`pf view-pocket`はPyMOL(`pymol-open-source`)を使う。rdkitのバージョン要件(`rdkit>=2026.03.5`)が
共通環境`pharmoforge`と競合するため、[`alignview`](../alignview/README.md)と同様に専用の
conda/mamba環境(既定名`pymol`)にインストールする前提。

```bash
mamba create -n pymol -c conda-forge pymol-open-source
```

`fpocket`(conda-forge)は`pf find-pocket`に必要で、`pharmoforge`環境に直接インストールする。

```bash
mamba install -n pharmoforge -c conda-forge fpocket
```

## 実装方針

- fpocketの実行(サブプロセス起動)・出力パース(`<stem>_info.txt`のポケット記述子、
  `pockets/pocket<N>_atm.pdb`のポケット周辺残基)は共通パッケージ[`src/pocket`](../src/pocket)で行う
  (アトミックな技術要素として独立)。
- fpocketは入力ファイルと同じディレクトリに出力を書き出すため、入力構造は`--outdir`にコピーして
  から実行する。ポケット周辺残基の抽出には常に`-w pdb`(PDB形式で出力)を指定し、
  [`structio`](../src/structio)(ProDy)で読み込む(fpocketが生成するmmCIFは`pdbx_PDB_model_num`列を
  欠くためProDyでパースできないことを確認済み)。
- `pockets.json`への整形・書き出し・読み込み、およびTSV(1行1ポケット)への整形・書き出しは
  本機能固有(`src/pocketfinder/report.py`。TSV化にはpandasの`DataFrame.to_csv(sep="\t")`を使う、
  [`scaffoldanalyzer`](../scaffoldanalyzer/README.md)と同様の方式)。
- ポケット可視化用のPyMOLスクリプト組み立ては`src/pocketfinder/view.py`。PyMOLの起動自体は共通
  パッケージ[`src/pymolrun`](../src/pymolrun)(`mamba run -n <env> pymol <script>`)を使う
  (元々`alignview`固有だったが、本機能でも同じ起動ロジックが必要になったためアトミックな技術要素
  として切り出した)。

## テスト

```bash
pytest tests/pocketfinder tests/pocket tests/pymolrun
```

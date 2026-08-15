# fetcher

ChEMBLやRCSB PDBなど外部データベースから、化合物・活性データや蛋白構造データを取得する機能。

## 使い方

```bash
pf fetch <データ種別>=<識別子> --output <出力ファイル>
```

### 活性データの取得(ChEMBL)

指定した標的蛋白について、pChEMBL値を持つ活性データを全件取得する。
化合物構造は[ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)に倣って標準化し(塩の除去等)、
標準化後の構造が同一の化合物についてはpChEMBL値をmean/median/sdに集約してTSVに書き出す。

```bash
pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
```

識別子には以下のいずれの形式を与えてもよい(自動判別して相互変換する)。

- UniProt entry name (例: `CDK4_HUMAN`)
- UniProt accession (例: `P11802`)
- ChEMBL target id (例: `CHEMBL331`)

出力TSVの列: `smiles`(標準化後のcanonical SMILES)、`_median`、`_mean`、`_sd`(集約対象が1件のみの場合は空)、`_n`(集約件数)。
行は`_median`の降順にソートされる。
SMILESがパースできない、またはpChEMBL値が欠損している記録は集約対象から除外する(除外数はログに表示)。

### 構造データの取得(RCSB PDB)

単体取得。`--type`(`cif`/`pdb`)でフォーマットを明示できる。省略時は出力ファイルの拡張子(`.cif` / `.pdb`)から推定する。

```bash
pf fetch structure=9CSK --type=cif --output data/9csk.cif
```

複数のPDB IDをカンマ区切りでまとめて取得することもできる(`structures=`、複数形)。
この場合`--output`は出力ディレクトリになり、各ファイルは`<出力ディレクトリ>/<PDB_ID>.<拡張子>`として保存される。
`--type`は必須(ディレクトリ名からはフォーマットを推定できないため)。

```bash
pf fetch structures=6P8F,7SJ3,9CSK --type pdb --output data
```

## 実装方針

- ChEMBLへのアクセスは `chembl_webresource_client` を使わず、`requests` による自前のAPI呼び出しで行う([`src/fetcher/chembl.py`](../src/fetcher/chembl.py))。
- 蛋白識別子の相互マッピング(UniProt entry name / accession / ChEMBL target id)は共通パッケージ [`src/idmap`](../src/idmap) で行う。
- 化合物構造の標準化は共通パッケージ [`src/molstd`](../src/molstd) で行う(fetcher専用ではなく横断的な技術要素として分離。`chembl_structure_pipeline` + RDKitを使用)。

## テスト

```bash
pytest tests/fetcher tests/idmap tests/molstd
```

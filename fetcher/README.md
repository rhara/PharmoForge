# fetcher

ChEMBLやRCSB PDBなど外部データベースから、化合物・活性データや蛋白構造データを取得する機能。

## 使い方

```bash
pf fetch <データ種別>=<識別子> --output <出力ファイル>
```

### 活性データの取得(ChEMBL)

指定した標的蛋白について、pChEMBL値を持つ活性データ(化合物構造・活性値)を全件取得しTSVに書き出す。

```bash
pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
```

識別子には以下のいずれの形式を与えてもよい(自動判別して相互変換する)。

- UniProt entry name (例: `CDK4_HUMAN`)
- UniProt accession (例: `P11802`)
- ChEMBL target id (例: `CHEMBL331`)

出力TSVの列: `molecule_chembl_id`, `canonical_smiles`, `target_chembl_id`, `target_pref_name`,
`standard_type`, `standard_relation`, `standard_value`, `standard_units`, `pchembl_value`,
`assay_chembl_id`, `assay_description`, `document_chembl_id`

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

## テスト

```bash
pytest tests/fetcher tests/idmap
```

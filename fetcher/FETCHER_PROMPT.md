# fetcher 実装記録

このドキュメントは `fetcher` 機能を再現するための仕様記録。

## 目的

外部データベース(ChEMBL、RCSB PDB等)から、化合物・活性データや蛋白構造データを取得する。
PharmoForgeにおける「データ収集」機能を担う。

## CLI仕様

```
pf fetch <データ種別>=<識別子> --output <出力ファイル>
```

- `<データ種別>` は現時点で `activities` (ChEMBL活性データ)、`structure` (RCSB PDB構造データ単体)、
  `structures` (RCSB PDB構造データ複数)をサポート。
- `--output` / `-o` は必須。出力ファイル名は都度明示的に指定する前提。
- `--type` (`cif`/`pdb`) は構造データのフォーマット指定に使う。

### activities: ChEMBL活性データ取得

```
pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
```

- `<識別子>` にはUniProt entry name / UniProt accession / ChEMBL target idのいずれを与えてもよい。
  `src/idmap` の `resolve_chembl_target_id()` が自動判別して相互変換する。
- ChEMBL REST API (`https://www.ebi.ac.uk/chembl/api/data/activity.json`) を `target_chembl_id` と
  `pchembl_value__isnull=false` で絞り込み、`page_meta.next` を辿って全件取得する。
  `chembl_webresource_client` パッケージは使わず、`requests` による自前実装とする。
- 出力はTSV。列は `src/fetcher/chembl.py` の `ACTIVITY_FIELDS` を参照
  (`molecule_chembl_id`, `canonical_smiles`, `standard_type`, `standard_value`, `pchembl_value` 等)。

### structure: RCSB PDB構造データ取得(単体)

```
pf fetch structure=9CSK --type=cif --output data/9csk.cif
```

- RCSB PDBのダウンロードエンドポイント `https://files.rcsb.org/download/{PDB_ID}.{fmt}` から取得する。
- フォーマット(`cif`/`pdb`)は `--type` で指定する。省略時は出力ファイルの拡張子から推定する(どちらもなければcif)。

### structures: RCSB PDB構造データ取得(複数)

```
pf fetch structures=6P8F,7SJ3,9CSK --type pdb --output data
```

- `<識別子>` はPDB IDをカンマ区切りで複数指定する。
- `--output` は出力ディレクトリになり、各ファイルは `<出力ディレクトリ>/<PDB_ID>.<fmt>` として保存される。
- `--type` は必須(出力先がディレクトリのためファイル拡張子からフォーマットを推定できない)。

## 実装ファイル

- `src/fetcher/chembl.py` — ChEMBL活性データ取得・TSV書き出し
- `src/fetcher/structures.py` — RCSB PDB構造データ取得
- `src/fetcher/cli.py` — `pf fetch` サブコマンド(データ種別のディスパッチ)
- `src/idmap/identifiers.py` — 蛋白識別子の相互マッピング(fetcher以外からも使う横断的な機能として分離)

## テスト

ネットワークアクセスは `unittest.mock` でモックし、実際の外部APIへは接続しない。

```bash
pytest tests/fetcher tests/idmap
```

## 動作例(サンプルデータ)

CDK4(ヒト、UniProt: P11802、ChEMBL: CHEMBL331)を題材にした活性データ取得例:

```bash
pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
```

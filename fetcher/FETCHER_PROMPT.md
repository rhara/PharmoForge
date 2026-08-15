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
  `structures` (RCSB PDB構造データ複数)、`structure-af` (AlphaFold DB構造データ単体)、
  `structures-af` (AlphaFold DB構造データ複数)をサポート。
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
  `chembl_webresource_client` パッケージは使わず、`requests` による自前実装とする
  (`src/chembl/activity.py` の `fetch_activities()`。fetcher固有ではなくアトミックな技術要素として独立)。
- 取得したactivity recordの `canonical_smiles` をRDKit標準の`rdMolStandardize`のみで標準化する
  (`src/molstd` の `standardize_smiles()`、`rdMolStandardize.Cleanup()` → `FragmentParent()` → `Uncharger().uncharge()`
  で塩・電荷を正規化。外部パッケージ`chembl_structure_pipeline`には非依存の自前実装)。
- 標準化後のSMILES(`smiles`)が同一の化合物ごとに `pchembl_value` を集約し、
  median/mean/sd(sample stdev、1件のみの場合は空)、集約件数(`_n`)を計算する
  (`src/fetcher/activities.py` の `standardize_and_aggregate()`)。
  SMILESがパースできない・pChEMBL値が欠損している記録は集約対象から除外する。
- 出力行は`_median`の降順にソートする。
- 出力はTSV。列は `src/fetcher/activities.py` の `AGGREGATED_FIELDS` を参照
  (`smiles`, `_median`, `_mean`, `_sd`, `_n`)。
  元の`molecule_chembl_id`・`target_chembl_id`・アッセイ関連の列(`standard_type`等)は出力しない。

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

### structure-af / structures-af: AlphaFold DB構造データ取得(単体/複数)

```
pf fetch structure-af=P61626 --type=cif --output data/P61626.cif
pf fetch structures-af=P61626,P11802 --type pdb --output data
```

- `<識別子>` はUniProt entry name(例: `TYK2_HUMAN`)またはaccession(例: `P29597`)を指定する
  (`src/idmap`の`resolve_uniprot_accession()`が自動判別してaccessionに解決する。`structures-af=`では
  カンマ区切りの各識別子を個別に解決する)。オプション体系(`--type`省略時の拡張子推定、複数形での
  出力ディレクトリ扱い等)はRCSB PDBの`structure`/`structures`と同じ。
- AlphaFold DBのpredictionエンドポイント(`https://alphafold.ebi.ac.uk/api/prediction/{accession}`)からダウンロードURL(`cifUrl`/`pdbUrl`)を都度解決する(バージョン番号をURLに固定しない)。

## 実装ファイル

- `src/fetcher/activities.py` — ChEMBL活性データの標準化+集約・TSV書き出し(fetcher固有の集計ロジック)
- `src/fetcher/cli.py` — `pf fetch` サブコマンド(データ種別のディスパッチ)
- `src/chembl/activity.py` — ChEMBL活性データの生取得(アトミックな技術要素として分離)
- `src/idmap/identifiers.py` — 蛋白識別子の相互マッピング(アトミックな技術要素として分離。`resolve_chembl_target_id()`は`activities=`、`resolve_uniprot_accession()`は`structure-af=`/`structures-af=`が使用)
- `src/molstd/standardize.py` — 化合物構造の標準化(アトミックな技術要素として分離。`rdMolStandardize`のみの自前実装)
- `src/rcsb/download.py` — RCSB PDB構造データ取得(アトミックな技術要素として分離)
- `src/afdb/download.py` — AlphaFold DB構造データ取得(アトミックな技術要素として分離)

## テスト

ネットワークアクセスは `unittest.mock` でモックし、実際の外部APIへは接続しない。

```bash
pytest tests/fetcher tests/idmap tests/chembl tests/molstd tests/rcsb tests/afdb
```

## 動作例(サンプルデータ)

CDK4(ヒト、UniProt: P11802、ChEMBL: CHEMBL331)を題材にした活性データ取得例:

```bash
pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
```

# fetcher 実装記録

このドキュメントは `fetcher` 機能を再現するための仕様記録。

## 目的

外部データベース(ChEMBL、RCSB PDB等)から、化合物・活性データや蛋白構造データを取得する。
PharmoForgeにおける「データ収集」機能を担う。

## CLI仕様

```
pf fetch <データ種別>=<識別子> --outdir <出力ディレクトリ>
```

- `<データ種別>` は現時点で `activity` (ChEMBL活性データ)、`structure` (構造データ、
  `,`区切りで複数可)をサポート。
- `--outdir` / `-o` は必須。出力先ディレクトリを指定する。ファイル名は識別子から自動的に決まる
  (都度明示的にファイル名を指定する必要はない)。
- `--type` (`cif`/`pdb`) は構造データのフォーマット指定に使う。省略時は`cif`。
- `--af` (フラグ、`structure=`専用): 取得元をRCSB PDBからAlphaFold DBに切り替える。
  `activity=`と併用するとエラー。

### activity: ChEMBL活性データ取得

```
pf fetch activity=CDK4_HUMAN --outdir data
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
- 出力ファイル名は `<出力ディレクトリ>/<識別子>_activity.tsv`(上記の例では`data/CDK4_HUMAN_activity.tsv`)。

### structure: 構造データ取得(単体・複数、RCSB PDB / AlphaFold DB)

```
pf fetch structure=9CSK --type=cif --outdir data
pf fetch structure=6P8F,7SJ3,9CSK --type pdb --outdir data
pf fetch structure=P61626 --af --type=cif --outdir data
pf fetch structure=CDK1_HUMAN --af --type=cif --outdir data/cdk1
pf fetch structure=P61626,CDK4_HUMAN --af --type pdb --outdir data
```

- `<識別子>` はカンマ区切りで複数指定できる(区切りがなければ単体扱い)。
- `--af`指定なし(既定): RCSB PDBのダウンロードエンドポイント
  `https://files.rcsb.org/download/{PDB_ID}.{fmt}` から取得する。`<識別子>`はPDB ID。
  出力ファイル名は `<出力ディレクトリ>/<識別子(大文字化)>.<fmt>`。
- `--af`指定あり: `<識別子>`にUniProt entry name(例: `CDK1_HUMAN`)またはaccession(例: `O14519`)を
  指定する。ダウンロード用accessionは`src/idmap`の`resolve_uniprot_accession()`が自動判別して解決する。
  AlphaFold DBのpredictionエンドポイント(`https://alphafold.ebi.ac.uk/api/prediction/{accession}`)から
  ダウンロードURL(`cifUrl`/`pdbUrl`)を都度解決する(バージョン番号をURLに固定しない)。
  出力ファイル名は識別子にentry nameを与えた場合でも常にUniProt accession(上記で解決したもの)
  に`_AF`を付けたものを使う: `<出力ディレクトリ>/<accession>_AF.<fmt>`
  (例: `pf fetch structure=CDK1_HUMAN --af --type=cif --outdir data/cdk1` →
  `data/cdk1/O14519_AF.cif`)。
- フォーマット(`cif`/`pdb`)は `--type` で指定する。省略時は`cif`。

## 実装ファイル

- `src/fetcher/activities.py` — ChEMBL活性データの標準化+集約・TSV書き出し(fetcher固有の集計ロジック)
- `src/fetcher/cli.py` — `pf fetch` サブコマンド(データ種別・`--af`フラグのディスパッチ)
- `src/chembl/activity.py` — ChEMBL活性データの生取得(アトミックな技術要素として分離)
- `src/idmap/identifiers.py` — 蛋白識別子の相互マッピング(アトミックな技術要素として分離。`resolve_chembl_target_id()`は`activity=`、`resolve_uniprot_accession()`は`structure=`の`--af`指定時が使用。ダウンロード用accessionの解決と、出力ファイル名(`<accession>_AF.<fmt>`)にも同じaccessionを使う)
- `src/molstd/standardize.py` — 化合物構造の標準化(アトミックな技術要素として分離。`rdMolStandardize`のみの自前実装)
- `src/rcsb/download.py` — RCSB PDB構造データ取得(アトミックな技術要素として分離)
- `src/afdb/download.py` — AlphaFold DB構造データ取得(アトミックな技術要素として分離)

## 設計変更の経緯

### `structure-af=`/`structures-af=` → `structure=`/`structures=` + `--af`

当初はRCSB PDBとAlphaFold DBを別々のデータ種別(`structure`/`structures` vs
`structure-af`/`structures-af`)として実装していたが、ユーザーからの要望により、
データ種別は`structure=`/`structures=`に統一し、取得元をRCSB PDB(既定)/AlphaFold DB
(`--af`フラグ指定時)で切り替える方式に変更した。オプション体系(`--type`、複数形での
出力ディレクトリ扱い等)は共通のまま、`--af`の有無で内部的に`rcsb.fetch_structure`と
`afdb.fetch_structure`のどちらを呼ぶかを切り替える(`--af`指定時のみ識別子を
`idmap.resolve_uniprot_accession()`で解決する)。`--af`は`activities=`と併用するとエラー。

### `activities=`/`structure=`+`structures=` + `--output` → `activity=`/`structure=`(単複統合) + `--outdir`

その後さらにユーザーからの要望により、以下のように変更した。

- `activities=`(複数形)→`activity=`に変更。1標的の活性データをまとめて取得する処理内容は
  変わらないが、データ種別名を単数形に統一した。
- `structure=`(単体)/`structures=`(複数)を`structure=`に統合。識別子を`,`区切りで複数
  指定すれば複数取得、区切りがなければ単体取得と、識別子の書き方だけで単複を判別する。
- `--output`/`-o`(出力ファイルパスを都度明示、複数形の場合のみディレクトリ扱い)を廃止し、
  `--outdir`/`-o`(常に出力先ディレクトリを指定)に統一した。ファイル名は識別子から自動的に
  決まる: `activity=`は`<識別子>_activity.tsv`、`structure=`は`<識別子>.<fmt>`
  (`--af`指定時はaccessionを与えても常にUniProt entry nameを使う。`idmap.resolve_uniprot_entry_name()`
  で解決)。
- `--type`の省略時デフォルトを、旧来の「出力ファイル拡張子からの推定(既定cif)」から
  単純な既定値`cif`に変更した(出力ファイル名を都度指定しなくなったため、拡張子からの
  推定という概念自体が不要になった)。

### `--af`指定時の出力ファイル名: UniProt entry name → UniProt accession + `_AF`

上記の変更直後は`--af`指定時の出力ファイル名にUniProt entry name(`idmap.resolve_uniprot_entry_name()`
で解決)を使っていたが、ユーザーからの指摘により、UniProt accessionに`_AF`を付けたもの
(`<accession>_AF.<fmt>`)に変更した。これに伴い`resolve_uniprot_entry_name()`/
`accession_to_entry_name()`は不要になったため`src/idmap`から削除し、`--af`時の識別子解決は
`resolve_uniprot_accession()`(ダウンロード用accessionの解決と出力ファイル名の両方に使う)のみに戻した。

### バグ修正: `entry_name_to_accession()`があいまい検索により誤ったaccessionを返すことがある問題

`idmap.identifiers.entry_name_to_accession()`は当初UniProtの検索エンドポイント
(`GET /uniprotkb/search?query=id:<entry_name>`)を使い、先頭の検索結果(`results[0]`)を
採用していたが、この検索はトークン化された全文検索でありentry nameの完全一致ではないため、
無関係のエントリが先頭に返ることがあった。実例: `id:CDK1_HUMAN`で検索すると、完全一致する
`CDK1_HUMAN`(P06493、Cyclin-dependent kinase 1)より先に`CDKA1_HUMAN`(O14519、CDK2-associated
protein 1、全く別の蛋白)が返り、`structure=CDK1_HUMAN --af`が誤ってO14519の構造をダウンロード
していた。これはユーザーからの指摘で発覚した重大なバグ。

修正: UniProtのエントリ取得エンドポイント(`GET /uniprotkb/{id}.json`、`src/uniprot/entry.py`の
`UNIPROT_ENTRY_URL`と同種)はentry name/accessionのどちらを与えても該当accessionへ一意に解決される
(存在しない識別子の場合はHTTP 400)ため、こちらを使うように変更した。あいまいさのある検索
エンドポイントは`idmap`からは使わないようにした。

## テスト

ネットワークアクセスは `unittest.mock` でモックし、実際の外部APIへは接続しない。

```bash
pytest tests/fetcher tests/idmap tests/chembl tests/molstd tests/rcsb tests/afdb
```

## 動作例(サンプルデータ)

CDK4(ヒト、UniProt: P11802、ChEMBL: CHEMBL331)を題材にした活性データ取得例:

```bash
pf fetch activity=CDK4_HUMAN --outdir data
```

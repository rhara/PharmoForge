# fetcher

ChEMBLやRCSB PDBなど外部データベースから、化合物・活性データや蛋白構造データを取得する機能。

## 使い方

```bash
pf fetch <データ種別>=<識別子> --outdir <出力ディレクトリ>
```

出力先は常に`--outdir`/`-o`で指定するディレクトリで、ファイル名は識別子から自動的に決まる。

### 活性データの取得(ChEMBL)

指定した標的蛋白について、pChEMBL値を持つ活性データを全件取得する。
化合物構造は共通パッケージ[`src/molstd`](../src/molstd)で標準化し(塩の除去等)、
標準化後の構造が同一の化合物についてはpChEMBL値をmean/median/sdに集約してTSVに書き出す。

```bash
pf fetch activity=CDK4_HUMAN --outdir data
```

上記は`data/CDK4_HUMAN_activity.tsv`(`<出力ディレクトリ>/<識別子>_activity.tsv`)として保存される。

識別子には以下のいずれの形式を与えてもよい(自動判別して相互変換する)。

- UniProt entry name (例: `CDK4_HUMAN`)
- UniProt accession (例: `P11802`)
- ChEMBL target id (例: `CHEMBL331`)

出力TSVの列: `smiles`(標準化後のcanonical SMILES)、`_median`、`_mean`、`_sd`(集約対象が1件のみの場合は空)、`_n`(集約件数)。
行は`_median`の降順にソートされる。
SMILESがパースできない、またはpChEMBL値が欠損している記録は集約対象から除外する(除外数はログに表示)。

### 構造データの取得(RCSB PDB / AlphaFold DB)

`--type`(`cif`/`pdb`)でフォーマットを指定する(省略時は`cif`)。
既定の取得元はRCSB PDB(識別子はPDB ID)。`--af`を付けるとAlphaFold DB(識別子はUniProt
entry name(例: `TYK2_HUMAN`)またはaccession(例: `P29597`)のいずれでもよい。ダウンロード用の
accessionは内部で`idmap.resolve_uniprot_accession`により解決される)から取得する。

識別子はカンマ区切りで複数指定できる(区切りがなければ単体扱い)。各ファイルは
`<出力ディレクトリ>/<識別子>.<拡張子>`として保存される。ただし`--af`指定時はファイル名に
UniProt accession(識別子にentry nameを与えた場合でも、内部で`idmap.resolve_uniprot_accession`
により解決したaccession)に`_AF`を付けたものを使う。

```bash
pf fetch structure=9CSK --type=cif --outdir data
pf fetch structure=6P8F,7SJ3,9CSK --type pdb --outdir data
pf fetch structure=P61626 --af --type=cif --outdir data
pf fetch structure=CDK1_HUMAN --af --type=cif --outdir data/cdk1
pf fetch structure=P61626,CDK4_HUMAN --af --type pdb --outdir data
```

`structure=P61626 --af`は`data/P61626_AF.cif`として、
`structure=CDK1_HUMAN --af`は`data/cdk1/O14519_AF.cif`(entry nameから解決したaccession)として保存される。

## 関数一覧

fetcher固有(アトミックでない、標準化後の値の集約・TSV書き出し等)の実装。
アトミックな共通パッケージ側の関数(ChEMBL活性データ生取得、RCSB構造ダウンロード等)は[API.md](../API.md)を参照。

### `src/fetcher/activities.py`

| 関数 | 説明 |
| --- | --- |
| `standardize_and_aggregate(records: list[dict]) -> list[dict]` | 取得した活性レコードの化合物構造を`molstd.standardize_smiles`で標準化し、同一構造ごとにpChEMBL値をmean/median/sdに集約する(`_median`降順)。 |
| `write_activities_tsv(records: list[dict], output: Path) -> None` | 集約済みレコードをTSV(`smiles, _median, _mean, _sd, _n`)として書き出す。 |

ChEMBL活性データの生取得(`fetch_activities`)、構造データ取得(`fetch_structure`/`fetch_structures`)は
fetcher固有ではなくアトミックな技術要素として、それぞれ[`src/chembl`](../src/chembl)、
[`src/rcsb`](../src/rcsb)に切り出し済み([API.md](../API.md)参照)。

## 実装方針

- ChEMBLへのアクセスは `chembl_webresource_client` を使わず、`requests` による自前のAPI呼び出しで行う(共通パッケージ [`src/chembl`](../src/chembl))。
- 蛋白識別子の相互マッピング(UniProt entry name / accession / ChEMBL target id)は共通パッケージ [`src/idmap`](../src/idmap) で行う。
- 化合物構造の標準化は共通パッケージ [`src/molstd`](../src/molstd) で行う(RDKit標準の`rdMolStandardize`のみを用いた自前実装、外部パッケージ`chembl_structure_pipeline`には非依存)。
- RCSB PDB構造ファイルのダウンロードは共通パッケージ [`src/rcsb`](../src/rcsb) で行う。
- AlphaFold DB構造ファイルのダウンロードは共通パッケージ [`src/afdb`](../src/afdb) で行う。

## テスト

```bash
pytest tests/fetcher tests/idmap tests/chembl tests/molstd tests/rcsb tests/afdb
```

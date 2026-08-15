# API

`src`以下のアトミックな(その機能固有のドメインロジックに閉じず、単体で再利用しうる)関数・パッケージの一覧。
記録対象は「現に複数の機能から使われているか」ではなく「再利用しうるアトミックな処理か」で判断する
(現時点の呼び出し元が1機能のみでも対象に含める、[CLAUDE.md](CLAUDE.md)参照)。
機能固有のドメインロジック(例: `src/fetcher`のChEMBL APIパラメータ組み立て)は各機能のREADME/PROMPT記録を参照。

## `src/core`

汎用ユーティリティ、および`pf`コマンド本体。

### `core.logging_utils`

| 関数 | 説明 |
| --- | --- |
| `get_logger(name: str) -> logging.Logger` | verboseな進捗表示用のロガーを返す。時刻・モジュール名付きのフォーマットで標準出力に流す。 |

### `core.cli`

`pf`の起点となるclickグループ`cli`。各機能パッケージが`cli.add_command(...)`でサブコマンドを登録する。

## `src/idmap`

蛋白識別子(UniProt entry name / UniProt accession / ChEMBL target id)の相互マッピング。

### `idmap.identifiers`

| 関数 | 説明 |
| --- | --- |
| `looks_like_chembl_id(identifier: str) -> bool` | ChEMBL target id(`CHEMBL\d+`)の形式かどうかを判定する。 |
| `looks_like_uniprot_accession(identifier: str) -> bool` | UniProt accessionの形式かどうかを判定する。 |
| `entry_name_to_accession(entry_name: str) -> str` | UniProt entry name(例: `CDK4_HUMAN`)をaccession(例: `P11802`)に変換する(UniProt REST APIを使用)。 |
| `accession_to_chembl_target_id(accession: str) -> str` | UniProt accessionをChEMBL target id(例: `CHEMBL331`)に変換する(ChEMBL REST APIを使用)。 |
| `resolve_chembl_target_id(identifier: str) -> str` | UniProt entry name / accession / ChEMBL target idのいずれを与えてもChEMBL target idを返す(上記3関数を組み合わせた入口)。 |

現時点では「識別子 → ChEMBL target id」の方向のみ実装。逆方向やUniProt accession単体の解決が他機能で必要になった時点で追加する。

## `src/chembl`

ChEMBL REST APIからの活性データ取得(`chembl_webresource_client`は使わず`requests`による自前実装)。

### `chembl.activity`

| 関数 | 説明 |
| --- | --- |
| `fetch_activities(target_chembl_id: str, page_size: int = 1000) -> list[dict]` | 指定したChEMBL target idについて、pChEMBL値を持つ活性データをChEMBL REST APIから全件取得する(ページネーション追従)。 |

## `src/molstd`

化合物構造の標準化。RDKit標準の`rdMolStandardize`のみを用いた自前実装
([ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)と同種の処理を、
外部パッケージ`chembl_structure_pipeline`には依存せず`Cleanup`/`FragmentParent`/`Uncharger`で組み立てている)。

### `molstd.standardize`

| 関数 | 説明 |
| --- | --- |
| `standardize_smiles(smiles: str) -> str \| None` | SMILESを`rdMolStandardize`(`Cleanup`→`FragmentParent`→`Uncharger`)で標準化し、親構造(塩等を除いた形)のcanonical SMILESを返す。パースに失敗した場合は`None`。 |

## `src/rcsb`

RCSB PDBからの構造ファイルダウンロード。

### `rcsb.download`

| 関数 | 説明 |
| --- | --- |
| `fetch_structure(pdb_id: str, output: Path, fmt: str \| None = None) -> None` | PDB IDの構造ファイルを1件ダウンロードし`output`に保存する。`fmt`(`cif`/`pdb`)省略時は`output`の拡張子から推定(既定`cif`)。 |
| `fetch_structures(pdb_ids: list[str], output_dir: Path, fmt: str) -> None` | 複数のPDB IDをまとめてダウンロードし、`output_dir/<PDB_ID>.<fmt>`として保存する。 |

## `src/afdb`

AlphaFold DBからの蛋白予測構造ファイルダウンロード。

### `afdb.download`

| 関数 | 説明 |
| --- | --- |
| `fetch_structure(accession: str, output: Path, fmt: str \| None = None) -> None` | UniProt accessionの構造ファイルを1件ダウンロードし`output`に保存する。`fmt`(`cif`/`pdb`)省略時は`output`の拡張子から推定(既定`cif`)。ダウンロードURLはAlphaFold DBのpredictionエンドポイント(`https://alphafold.ebi.ac.uk/api/prediction/{accession}`)から都度解決する(バージョンをURLに固定しない)。 |
| `fetch_structures(accessions: list[str], output_dir: Path, fmt: str) -> None` | 複数のUniProt accessionをまとめてダウンロードし、`output_dir/<ACCESSION>.<fmt>`として保存する。 |

## `src/molscaffold`

化合物のBemis-Murckoスキャフォールド計算。

### `molscaffold.scaffold`

| 関数 | 説明 |
| --- | --- |
| `compute_scaffold(smiles: str) -> str \| None` | SMILESからBemis-Murckoスキャフォールド(canonical SMILES)を求める(RDKitの`MurckoScaffold.GetScaffoldForMol`)。パースに失敗した場合は`None`。 |

## `src/actbin`

活性値(または任意の連続値)の分位点ビニング。

### `actbin.binning`

| 関数 | 説明 |
| --- | --- |
| `assign_activity_bins(df: pd.DataFrame, activity_col: str, high_quantile: float = 0.75, low_quantile: float = 0.25) -> pd.DataFrame` | `activity_col`の分位点で各行を`high`/`mid`/`low`に分類した`bin`列を追加する。`low_quantile < high_quantile`(0〜1)でない場合は`ValueError`。 |

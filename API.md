# API

`src`以下の汎用的な(複数の機能から使われうる)関数・パッケージの一覧。
機能固有の実装(例: `src/fetcher`)は各機能のREADME/PROMPT記録を参照。

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

## `src/molstd`

化合物構造の標準化([ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)に倣う)。

### `molstd.standardize`

| 関数 | 説明 |
| --- | --- |
| `standardize_smiles(smiles: str) -> str \| None` | SMILESを`chembl_structure_pipeline`で標準化し、親構造(塩等を除いた形)のcanonical SMILESを返す。パースに失敗した場合は`None`。 |

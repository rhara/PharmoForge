# API

`src`以下の、再利用しうる関数・パッケージの一覧。ノートブック(`notebooks/`)はこれらの関数の実例集
という位置づけであり、関数のシグネチャ・挙動の正はここに置く。
記録対象は「現に複数の機能から使われているか」ではなく「再利用しうる処理か」で判断する
(現時点の呼び出し元が1機能・1ノートブックのみでも対象に含める、[CLAUDE.md](CLAUDE.md)参照)。
root直下にREADME/PROMPT記録を持つ機能パッケージ(`proteinprep`/`docking`等)であっても、含まれる関数が
再利用しうるならAPI.mdにも記録する(root README/PROMPTの有無はAPI.md記載の可否を決めない。両者は
役割が異なり、README/PROMPTは機能としての使い方・経緯・環境構築、API.mdは関数シグネチャ・挙動の
リファレンス)。機能固有のドメインロジックそのもの(例: `src/fetcher`のChEMBL APIパラメータ組み立て)は
関数として切り出しにくいためAPI.md対象外とし、各機能のREADME/PROMPT記録側に置く。

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
| `resolve_uniprot_accession(identifier: str) -> str` | UniProt entry name / accessionのいずれを与えてもUniProt accessionを返す(`fetcher`の`structure=` `--af`指定時、ダウンロード用accessionと出力ファイル名(`<accession>_AF.<fmt>`)の両方の解決に使用)。 |
| `pdb_id_to_uniprot_accessions(pdb_id: str) -> list[str]` | PDB ID配下の全ポリマーエンティティに紐づくUniProt accessionを重複なく出現順で返す(RCSB Data API `data.rcsb.org`を使用)。複合体で複数の蛋白質が含まれる場合は複数返る。`fetcher`の`structure=` `--type=fasta`(`--af`なし)時、PDB IDからUniProt accessionを解決するのに使用。 |

## `src/chembl`

ChEMBL REST APIからの活性データ取得(`chembl_webresource_client`は使わず`requests`による自前実装)。

### `chembl.activity`

| 関数 | 説明 |
| --- | --- |
| `fetch_activities(target_chembl_id: str, page_size: int = 1000) -> list[dict]` | 指定したChEMBL target idについて、pChEMBL値を持つ活性データをChEMBL REST APIから全件取得する(ページネーション追従)。 |

### `chembl.local`

ChEMBL Web APIが障害・レート制限等で使えない場合のフォールバック。ChEMBL公式配布のSQLiteデータベース
(`chembl_XX.db`)に対し`sqlite3`で直接クエリする。

| 関数 | 説明 |
| --- | --- |
| `resolve_target_chembl_id(accession: str, db_path: str \| Path) -> str` | UniProt accessionから、SINGLE PROTEINターゲットのChEMBL target idを解決する(`idmap.accession_to_chembl_target_id`のWeb API版と同じ入出力)。見つからなければ`ValueError`。 |
| `fetch_activities(target_chembl_id: str, db_path: str \| Path) -> list[dict]` | 指定したChEMBL target idについて、pChEMBL値を持つ活性データを取得する。各要素は`chembl.activity.fetch_activities`(Web API版)と同じ主要フィールド(`molecule_chembl_id`/`molecule_pref_name`/`canonical_smiles`/`standard_type`/`standard_value`/`standard_units`/`pchembl_value`/`assay_chembl_id`/`document_chembl_id`)を持つ。 |

### `chembl.aggregate`

複数のChEMBL標的にまたがる活性データの集計(化合物×標的単位の要約、化合物単位のロールアップ、
潜在活性化合物の抽出)。標的ごとの活性データ取得そのものは`chembl.local.fetch_activities`に委譲する。

| 関数 | 説明 |
| --- | --- |
| `collect_standardized_activities(targets_df: pd.DataFrame, db_path: str \| Path, target_id_col: str = "chembl_target_id", accession_col: str = "accession", entry_name_col: str = "entry_name") -> pd.DataFrame` | 複数のChEMBL標的について活性データを収集し、標準化SMILES・pChEMBL値のみのフラットなテーブル(列: `smiles`/`accession`/`entry_name`/`pchembl_value`)にする。`target_id_col`がNaNの行はスキップ。pChEMBL値・SMILESを欠く活性、標準化(`molstd.standardize_smiles`)に失敗した化合物は除外する。 |
| `summarize_compound_target_activity(activity_df: pd.DataFrame, group_cols: tuple[str, ...] = ("smiles", "accession", "entry_name"), value_col: str = "pchembl_value") -> pd.DataFrame` | 化合物×標的の単位で活性値のmedian/mean/std/個数を集計する(median降順)。 |
| `rollup_compound_summary(activity_summary_df: pd.DataFrame, smiles_col: str = "smiles", entry_name_col: str = "entry_name", value_col: str = "median") -> pd.DataFrame` | `summarize_compound_target_activity`の出力を化合物単位にロールアップする(`target_count`・`best_target_<entry_name_col>`・`best_<value_col>`、`best_<value_col>`降順)。 |
| `select_high_potency_compounds(df: pd.DataFrame, potency_col: str, potency_cutoff: float, mol_weight_col: str \| None = None, mol_weight_range: tuple[float, float] \| None = None) -> pd.DataFrame` | `potency_col >= potency_cutoff`の行を抽出する。`mol_weight_col`/`mol_weight_range`を両方指定した場合は分子量(drug-like範囲等)でも絞り込む。 |

## `src/molstd`

化合物構造の標準化。RDKit標準の`rdMolStandardize`のみを用いた自前実装
([ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)と同種の処理を、
外部パッケージ`chembl_structure_pipeline`には依存せず`Cleanup`/`FragmentParent`/`Uncharger`で組み立てている)。

### `molstd.standardize`

| 関数 | 説明 |
| --- | --- |
| `standardize_smiles(smiles: str) -> str \| None` | SMILESを`rdMolStandardize`(`Cleanup`→`FragmentParent`→`Uncharger`)で標準化し、親構造(塩等を除いた形)のcanonical SMILESを返す。パースに失敗した場合は`None`。 |

### `molstd.descriptors`

| 関数 | 説明 |
| --- | --- |
| `calc_mol_weight(smiles: str) -> float \| None` | SMILESから平均分子量(Da、RDKitの`Descriptors.MolWt`)を計算する。パースに失敗した場合は`None`。 |

## `src/rcsb`

RCSB PDBからの構造ファイルダウンロード・メタデータ取得。

### `rcsb.download`

| 関数 | 説明 |
| --- | --- |
| `fetch_structure(pdb_id: str, output: Path, fmt: str \| None = None) -> None` | PDB IDの構造ファイルを1件ダウンロードし`output`に保存する。`fmt`(`cif`/`pdb`)省略時は`output`の拡張子から推定(既定`cif`)。 |
| `fetch_structures(pdb_ids: list[str], output_dir: Path, fmt: str) -> None` | 複数のPDB IDをまとめてダウンロードし、`output_dir/<PDB_ID>.<fmt>`として保存する。 |

### `rcsb.metadata`

| 関数 | 説明 |
| --- | --- |
| `fetch_entry_info(pdb_id: str) -> dict` | PDB IDの実験手法・解像度等をRCSB Data APIから取得する。`{"pdb_id", "method", "resolution"}`を返す。`resolution`はX線構造以外では`None`になりうる。`method`の値はRCSBの表記(例: `"X-ray"`、`"EM"`)をそのまま返す。 |
| `fetch_entries_info(pdb_ids: list[str]) -> list[dict]` | 複数PDB IDのメタデータをまとめて取得する。 |

## `src/afdb`

AlphaFold DBからの蛋白予測構造ファイルダウンロード。

### `afdb.download`

| 関数 | 説明 |
| --- | --- |
| `fetch_structure(accession: str, output: Path, fmt: str \| None = None) -> None` | UniProt accessionの構造ファイルを1件ダウンロードし`output`に保存する。`fmt`(`cif`/`pdb`/`fasta`)省略時は`output`の拡張子から推定(既定`cif`)。`cif`/`pdb`のダウンロードURLはAlphaFold DBのpredictionエンドポイント(`https://alphafold.ebi.ac.uk/api/prediction/{accession}`)から都度解決する(バージョンをURLに固定しない)。`fasta`指定時はAlphaFold DBへは問い合わせず、[`src/uniprot`](#srcuniprot)の`fetch_fasta()`でUniProt本体から正規配列を直接取得する(AlphaFold DBのモデルは断片のことがあるため)。 |
| `fetch_structures(accessions: list[str], output_dir: Path, fmt: str) -> None` | 複数のUniProt accessionをまとめてダウンロードし、`output_dir/<ACCESSION>.<fmt>`として保存する。 |

## `src/uniprot`

UniProtエントリの取得と、創薬(構造生物学・メディシナルケミストリー)向け情報の抽出。

### `uniprot.entry`

| 関数 | 説明 |
| --- | --- |
| `fetch_entry(accession: str) -> dict` | UniProt accessionの生エントリJSON(UniProt REST API)を取得する。 |
| `fetch_entry_names(accessions: list[str]) -> dict[str, str]` | 複数のUniProt accessionのentry name(例: `P11802` -> `CDK4_HUMAN`)をバッチ取得エンドポイント(`/uniprotkb/accessions`、最大100件/リクエスト)でまとめて取得する。accessionごとに`fetch_entry`を呼ぶより大幅に少ないリクエスト数で済む。見つからなかったaccessionは返り値に含まれない。 |
| `fetch_fasta(accession: str) -> bytes` | UniProt accessionの正規配列をUniProt標準のFASTA形式(`https://rest.uniprot.org/uniprotkb/{accession}.fasta`)で取得する。 |
| `extract_protein_info(entry: dict) -> dict` | 生エントリJSONから創薬関連情報を平坦なdictに整理する。抽出項目: `accession`/`entry_name`/`protein_name`/`gene_name`/`organism`/`taxon_id`/`sequence`/`length`/`mol_weight`/`ec_numbers`/`function`/`keywords`/`diseases`/`active_sites`/`binding_sites`/`disulfide_bonds`/`glycosylation_sites`/`modified_residues`/`transmembrane_regions`/`signal_peptide`/`domains`/`pdb_structures`(id/method/resolution)/`alphafold_id`。 |
| `fetch_protein_info(accession: str) -> dict` | `fetch_entry` + `extract_protein_info` をまとめた入口。 |
| `parse_resolution(resolution: str \| None) -> float \| None` | `pdb_structures`の解像度表記(例: `"2.25 A"`)をÅ単位の`float`に変換する。NMR構造等で解像度が無い場合や数値として解釈できない場合は`None`。 |

## `src/blastsearch`

NCBI BLAST Web API(QBLAST)を用いた配列相同性検索。

### `blastsearch.ncbi`

| 関数 | 説明 |
| --- | --- |
| `submit_blast(sequence: str, program: str = "blastp", database: str = "pdb", entrez_query: str \| None = None) -> str` | BLASTジョブを投函し、Request ID(RID)を返す。`database`は`"pdb"`(PDBエントリの配列)、`"swissprot"`(UniProtKB/Swiss-Prot)等、NCBI BLASTが受け付ける値を指定できる。`entrez_query`はNCBI Entrezのクエリ構文(例: `"Homo sapiens[Organism]"`)で検索対象を絞り込む。 |
| `wait_for_blast(rid: str, poll_interval: float = 10.0, timeout: float = 600.0) -> None` | BLASTジョブの完了(`Status=READY`)までポーリングして待機する。失敗ステータスで`RuntimeError`、タイムアウトで`TimeoutError`。 |
| `fetch_hits(rid: str) -> list[dict]` | 完了したBLASTジョブのヒット一覧を取得する。各要素は`{"subject_id", "identity", "align_length", "evalue", "bit_score"}`。`subject_id`はデータベースに応じた生の識別子文字列で、`parse_pdb_subject_id`/`parse_uniprot_subject_id`で解釈する。 |
| `blast_search(sequence: str, program: str = "blastp", database: str = "pdb", entrez_query: str \| None = None, poll_interval: float = 10.0, timeout: float = 600.0) -> list[dict]` | `submit_blast` + `wait_for_blast` + `fetch_hits` をまとめた入口。 |
| `parse_pdb_subject_id(subject_id: str) -> tuple[str, str]` | `database="pdb"`のヒットのsubject id(例: `pdb\|6GZM\|A`、`6GZM_A`)からPDB IDとchainを取り出す。 |
| `parse_uniprot_subject_id(subject_id: str) -> str` | `database="swissprot"`等のヒットのsubject id(例: `sp\|P11802\|CDK4_HUMAN`、`P11802.1`)からUniProt accessionを取り出す(バージョン番号は除く)。 |
| `best_hit_per_accession(hits: list[dict], exclude_accession: str \| None = None) -> list[dict]` | UniProt accessionごとに最良ヒット(evalue最小)だけを残し、evalue昇順で返す(各ヒットに`"accession"`キーを追加)。`exclude_accession`(クエリ自身のaccession等)は除外できる。 |
| `format_evalue(evalue: float) -> str` | e-valueを表示用に整形する。0はそのまま`"0"`、それ以外は有効数字2桁の指数表記(例: `9.47e-95` -> `"9.5e-95"`)。 |

### `blastsearch.cache`

BLASTジョブのファイルキャッシュ・再開可能な実行(notebookでの繰り返し実行を想定)。

| 関数 | 説明 |
| --- | --- |
| `run_cached_blast(sequence: str, cache_dir: Path, program: str = "blastp", database: str = "pdb", entrez_query: str \| None = None, poll_interval: float = 10.0, timeout: float = 600.0) -> list[dict]` | `blast_search`にファイルキャッシュ(`cache_dir/blast_hits.pkl`)と再開可能性を加える。既にRID(`cache_dir/blast_rid.txt`)が投函済みならジョブを再投函せず待機を再開する。待機が`TimeoutError`の場合はRIDキャッシュを残したまま伝播(再度呼べば再開できる)、ジョブ失敗(`RuntimeError`)の場合はRIDキャッシュを破棄する(再投函が必要なため)。 |

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

## `src/structio`

PDB/CIF構造ファイルの読み書き(拡張子で自動判別、[ProDy](http://prody.csb.pitt.edu/)の
`AtomGroup`を介する)。`structfit`・`proteinextract`・`alignview`・`sequencealign`が共通で利用する。

### `structio.io`

| 関数 | 説明 |
| --- | --- |
| `parse_structure(path: Path) -> Atomic` | PDB/CIF形式の構造ファイルを拡張子(`.pdb`/`.cif`/`.mmcif`)で自動判別して読み込む。CIFは`unite_chains=True`で読み込み、チェーンIDに`auth_asym_id`(PyMOLやRCSB Webサイトで見えるチェーンID)を使う(ProDyの既定`label_asym_id`は同じauthチェーンに属する水分子・リガンド等を別チェーンに細分化するため)。 |
| `write_structure(atoms: Atomic, path: Path) -> None` | ProDyの`AtomGroup`(または`select()`の結果)を、拡張子で自動判別したPDB/CIF形式で書き出す。出力先ディレクトリが存在しない場合は作成する。 |

### `structio.resolve`

複数構造ファイルをコマンドライン引数から解決する共通ロジック(`--indir`・拡張子省略への対応)。
`alignview`・`sequencealign`のCLIが共用する。

| 関数 | 説明 |
| --- | --- |
| `resolve_structure_tokens(tokens: tuple[str, ...], extensions: tuple[str, ...] = (".cif", ".mmcif", ".pdb")) -> list[Path]` | `click`の`UNPROCESSED`引数列を構造ファイルパスのリストに解決する。`--indir DIR`は繰り返し指定可能で、以降のファイル名(拡張子省略可、`extensions`の順で試す)を`DIR`配下から解決する。`extensions`は呼び出し元ごとに変更でき、`sequencealign`は`(".cif", ".mmcif", ".pdb", ".fasta")`を渡す(`alignview`はPyMOLでの3次元表示が前提のため既定のまま)。`/`を含む指定(または絶対パス)は`--indir`によらずカレントディレクトリ相対 or 絶対パスとして扱う。トークンが空、または`--indir`に値がない場合は`click.UsageError`。 |

### `structio.fasta`

| 関数 | 説明 |
| --- | --- |
| `parse_fasta(path: Path) -> list[tuple[str, str]]` | FASTAファイルをBiopython `SeqIO`で読み込み、`(ヘッダー, 配列)`のタプルのリストを返す。`sequencealign`が`.fasta`入力の配列読み込みに使う。 |

## `src/structfit`

同一蛋白の複数構造間で、残基番号の対応のみに基づく剛体重ね合わせ(rigid-body fit)を計算する
([ProDy](http://prody.csb.pitt.edu/)を使用、conda-forgeでのビルド提供が不安定なためpipでインストール)。
配列アラインメントを一切行わないため、構造間でPDBの残基番号(UniProt基準等)が揃っていることが前提。
構造の読み込みは[`structio`](#srcstructio)を利用する。

### `structfit.fit`

| 関数 | 説明 |
| --- | --- |
| `fit_by_residue_number(mobile_path: Path, target_path: Path) -> FitResult` | `mobile_path`・`target_path`の構造(PDB/CIF)を読み込み、残基番号が一致するCA原子同士を対応付けて`mobile`を`target`に重ね合わせる剛体変換を求める。両構造とも複数鎖を含む場合、共通残基番号数が最大になる鎖の組を自動選択する(NCS等による複数コピー対策)。共通の残基番号(CA)が1つもない場合は`ValueError`。 |
| `fit_by_residue_pairs(mobile_path: Path, target_path: Path, mobile_chain: str, target_chain: str, resnum_pairs: list[tuple[int, int]]) -> FitResult` | `fit_by_residue_number`の、残基番号体系が異なる蛋白間版。呼び出し側が用意した`resnum_pairs`((mobile側残基番号, target側残基番号)の対応リスト、例: 配列アラインメントで対応付けたポケット周辺残基)のCA原子同士を対応付けて重ね合わせる。対応するCA原子が1組もない場合は`ValueError`。 |
| `apply_fit(fit_result: FitResult, atoms: Atomic) -> Atomic` | `fit_result`の剛体変換を`atoms`(そのmobile構造由来のAtomic、選択結果でも可)に適用する。選択結果を渡しても親AtomGroup全体の座標が変換される(チェーンID等のラベルは変更しない)。 |
| `FitResult` | `fit_by_residue_number`/`fit_by_residue_pairs`の戻り値(dataclass)。`matrix`(4x4 numpy配列、行優先。`v' = matrix @ [x, y, z, 1]`)、`rmsd: float`、`n_residues: int`、`mobile_chain: str`、`target_chain: str`。 |

### `structfit.chainselect`

| 関数 | 説明 |
| --- | --- |
| `find_best_chain_for_residues(mobile_chains: list[ChainSequence], reference_sequence: str, reference_resnums: list[int]) -> BestChainCoverage \| None` | `mobile_chains`(`seqextract.get_chain_sequences`の結果)の各チェーンを`reference_sequence`にアラインメント(`seqalign.align_to_reference`)し、`reference_resnums`(基準配列側の残基番号集合、例: ポケット周辺残基)を最も多くカバーするチェーンを選ぶ。結晶構造のNCSコピーや無関係な鎖が混在する場合でも、目的の残基集合を最もよくカバーするチェーンを自動選択できる。1件もカバーしない場合は`None`。 |
| `BestChainCoverage` | 選定結果(dataclass)。`chain_id: str`、`resnum_pairs: list[tuple[int, int]]`(`fit_by_residue_pairs`にそのまま渡せる形式)。 |

## `src/seqextract`

構造(ProDy Atomic)から蛋白チェーンごとの配列(1文字表記)+残基番号を抽出する。CA原子(観測された
残基のみ)に基づく配列であり、電子密度が見えず欠損した残基は含まれない(UniProtの完全配列とは
異なりうる)。`sequencealign`向けにアトミックな技術要素として切り出した。

### `seqextract.chains`

| 関数 | 説明 |
| --- | --- |
| `get_chain_sequences(atoms: Atomic) -> list[ChainSequence]` | `atoms`に含まれる蛋白チェーンごとに配列を抽出する(チェーンID昇順)。水分子等、蛋白でないチェーン(CA原子を持たないもの)は結果に含めない。 |
| `ChainSequence` | 1蛋白チェーン分の配列情報(dataclass)。`chain_id: str`、`sequence: str`(1文字表記)、`resnums: list[int]`(`sequence[i]`に対応する残基番号)、`length`(プロパティ、配列長)。 |

## `src/structcompare`

構造(ProDy Atomic)間のチェーン単位配列比較。ProDyの`matchChains`をラップする。`matchChains`は
まず残基番号・残基名による直接対応付け(`pwalign=False`)を試み、失敗時にBiopythonによる配列
アラインメントにフォールバックする(`pwalign=True`)。現行のprodyバージョン(2.6.1)では
`pwalign=True`時に返る`AtomMap`の残基番号が実際の対応関係を反映しないことを実データで確認して
いるため(重ね合わせ後RMSDが数十Åに達する)、残基単位の対応付けには`pwalign=False`のみを使う
(詳細は[`sequencealign/SEQUENCEALIGN_PROMPT.md`](sequencealign/SEQUENCEALIGN_PROMPT.md)参照)。
`sequencealign`向けにアトミックな技術要素として切り出した。

### `structcompare.compare`

| 関数 | 説明 |
| --- | --- |
| `match_chains(atoms_a: Atomic, atoms_b: Atomic) -> list[ChainMatch]` | `atoms_a`と`atoms_b`のチェーンを全対全で比較し、対応が取れた組を%identity降順で返す(`pwalign=True`で広く候補を探す。参考情報用途、残基単位の対応には使わない)。 |
| `find_substitutions(atoms_a: Atomic, atoms_b: Atomic) -> SubstitutionReport` | `atoms_a`と`atoms_b`の最良一致チェーン同士で、残基番号ベース(`pwalign=False`)の対応付けが取れた範囲のアミノ酸置換を列挙する。対応が取れない場合は`matched=False`。 |
| `ChainMatch` | `chain_id_a: str`、`chain_id_b: str`、`seqid: float`(%)、`overlap: float`(%)、`n_matched: int`。 |
| `ResidueSubstitution` | `resnum: int`(a側・b側で共通の残基番号)、`resname_a: str`、`resname_b: str`。 |
| `SubstitutionReport` | `matched: bool`、`chain_id_b: str \| None`(atoms_b側で対応が取れたチェーンID)、`seqid: float \| None`、`overlap: float \| None`、`substitutions: list[ResidueSubstitution]`。欠損領域(a側にのみ存在する残基番号)は呼び出し側(`sequencealign`)が`chain_id_b`を使って`seqextract`の残基番号集合の差分から算出する。 |

## `src/seqalign`

構造を伴わない任意のアミノ酸配列(1文字表記)同士のペアワイズグローバルアラインメント。比較対象の
一方にProDyの`Atomic`が存在しない場合(例: UniProt正規配列やユーザー指定の基準配列との比較)に、
Biopythonの`Bio.Align.PairwiseAligner`を直接用いる(`structcompare`とは異なりProDyを介さない)。
`sequencealign`向けにアトミックな技術要素として切り出した。

### `seqalign.pairwise`

| 関数 | 説明 |
| --- | --- |
| `align_to_reference(ref_sequence: str, query_sequence: str, query_resnums: list[int]) -> AlignmentResult` | `ref_sequence`と`query_sequence`(長さは`query_resnums`と同じ)をBLOSUM62によるグローバルアラインメントし、%identity・%coverage・アミノ酸置換・欠失/挿入領域(`gaps`)を返す。 |
| `SequenceSubstitution` | `ref_pos: int`(基準配列内の位置、1始まり)、`ref_aa: str`、`query_resnum: int`(比較対象側の実際の残基番号)、`query_aa: str`。 |
| `SequenceGap` | 基準配列に対する欠失(`kind="deletion"`)またはqueryにのみ存在する挿入(`kind="insertion"`)領域。`ref_start`/`ref_end: int \| None`(基準配列内の範囲、insertionの場合None)、`length: int`、`before_query_resnum`/`after_query_resnum: int \| None`(query側でこの領域の直前・直後にある残基番号。配列末端の場合None)。 |
| `AlignmentResult` | `identity: float`(%、アラインメントされた位置のうち一致した割合)、`coverage: float`(%、`ref_sequence`全体のうちアラインメントされた位置の割合)、`aligned_length: int`、`substitutions: list[SequenceSubstitution]`、`gaps: list[SequenceGap]`、`query_by_ref_pos: dict[int, str]`(基準配列内の位置(1始まり)→対応するquery側のアミノ酸。ギャップの位置は含まない。`sequencealign`が配列アラインメントベースの整列表示・identityグリッド構築に使用)。 |

## `src/pocket`

[fpocket](https://github.com/Discngine/fpocket)(conda-forge、`pharmoforge`環境に別途インストールが
必要)の実行と出力パース。ポケット検出そのものはfpocketに委譲し、`<stem>_info.txt`(ポケット記述子)・
`pockets/pocket<N>_atm.pdb`(ポケットに面する原子)の読み取りのみを行う。`pocketfinder`向けに
アトミックな技術要素として切り出した。構造読み込みは[`structio`](#srcstructio)を利用する。

### `pocket.fpocket`

| 関数 | 説明 |
| --- | --- |
| `run_fpocket(structure_path: Path, work_dir: Path) -> list[Pocket]` | `structure_path`(PDB/CIF)を`work_dir`にコピーして`fpocket -f <コピー> -w pdb`を実行し、検出されたポケットをスコア降順で返す。fpocketは入力と同じディレクトリに`<stem>_out/`を書き出すため、`work_dir`にはfpocketの生出力(PyMOL/VMD可視化スクリプト等)がそのまま残る。ポケット周辺残基の抽出には常に`-w pdb`を指定する(fpocket生成のmmCIFは`pdbx_PDB_model_num`列を欠きProDyでパース不可のため)。fpocketが非0で終了した場合は`RuntimeError`。 |
| `Pocket` | 1ポケット分の情報(dataclass)。`pocket_id: int`、`score: float`、`druggability_score: float`、`n_alpha_spheres: int`、`volume: float`、`total_sasa`/`polar_sasa`/`apolar_sasa: float`、`hydrophobicity_score: float`、`residues: list[PocketResidue]`。 |
| `PocketResidue` | ポケットに面する残基1件(dataclass)。`chain_id: str`、`resnum: int`(auth番号)、`resname: str`。 |

### `pocket.selection`

| 関数 | 説明 |
| --- | --- |
| `select_pocket_by_anchor_overlap(pockets: list[Pocket], anchor_resnums: set[int], chain_id: str) -> PocketSelection` | `anchor_resnums`(保存モチーフ・相同蛋白のリガンド接触残基等、生物学的根拠のある残基集合)との重なりが最大のポケットを選ぶ。fpocketのスコア(druggability等)は「目的のポケットかどうか」を直接表さないため、スコアではなく既知のアンカー残基をどれだけ含むかで選ぶ。重なりが0件のポケットしかない場合は`ValueError`。 |
| `PocketSelection` | 選定結果(dataclass)。`pocket: Pocket`、`overlap_resnums: list[int]`、`overlap: int`(プロパティ、`len(overlap_resnums)`)。 |

## `src/kinasemotifs`

プロテインキナーゼドメインの保存モチーフ(P-loop/触媒Lys/HRD/DFG)を配列から検出する。ATP結合部位を
構成することが構造生物学的に確立している位置をアンカーとして使う用途(例: fpocket検出ポケットの
中からATP結合部位を特定する)を想定。CDK20の調査で`docking`向けにATP結合部位を特定する過程で、
どのキナーゼにも使える汎用ロジックとして切り出した。

### `kinasemotifs.motifs`

| 関数 | 説明 |
| --- | --- |
| `find_kinase_motifs(sequence: str) -> KinaseMotifs` | タンパク質配列からP-loop(`G.G..G`)・触媒Lys(`VA[LIV]K`)・HRD(`HRD`)・DFG(`DFG`)・DFG+1(DFGの直後、多くのキナーゼでback pocketの壁を構成することが知られる位置)を検出する。全て見つからない場合は`ValueError`。個々のモチーフが検出できないこと自体は許容する(該当フィールドが`None`になる)。 |
| `KinaseMotifs` | 検出結果(dataclass、1始まりの残基番号、区間は両端含む)。`p_loop`/`hrd`/`dfg: tuple[int, int] \| None`、`catalytic_lys`/`dfg_plus1: int \| None`、`anchor_resnums: set[int]`(プロパティ、検出できた全モチーフの残基番号をまとめた集合)。 |

## `src/ligandcontacts`

相同蛋白の複数のX線構造にまたがる、共結晶化リガンドの接触残基のコンセンサスを求める。単一構造だけの
接触では採用せず、一定割合以上の構造で再現された残基だけを採用することで、フラグメントスクリーニング
のオフターゲットヒットや無関係な部位への結合の影響を抑える。CDK20の調査でATP結合部位を相同蛋白の
共結晶化リガンドから特定する過程で、蛋白・部位によらず使える汎用ロジックとして切り出した。
構造読み込みは[`structio`](#srcstructio)、鎖配列抽出は[`seqextract`](#srcseqextract)、基準配列への
番号マッピングは[`seqalign`](#srcseqalign)を利用する。

### `ligandcontacts.consensus`

| 関数 | 説明 |
| --- | --- |
| `find_consensus_ligand_contacts(structure_paths: list[Path], reference_sequence: str, contact_distance: float = 4.5, min_ligand_atoms: int = 8, min_fraction: float = 0.2, excluded_resnames: frozenset[str] = DEFAULT_EXCLUDED_RESNAMES) -> ConsensusLigandContacts` | `structure_paths`(相同蛋白のPDB/CIF)に含まれる共結晶化リガンドの接触残基(CA、`contact_distance`Å以内)を、配列アラインメントで`reference_sequence`の番号にマッピングし、`min_fraction`(リガンドを含む構造数に対する割合、最低2構造は要求)以上で再現された残基だけをコンセンサスとして返す。`min_ligand_atoms`未満のヘテロ残基・`excluded_resnames`(結晶化添加物・修飾残基等)は無視する。 |
| `ConsensusLigandContacts` | 集計結果(dataclass)。`anchor_resnums: list[int]`(コンセンサス残基)、`n_ligands: int`、`contact_counts: dict[int, int]`(閾値適用前の生の接触回数)、`min_count: int`(採用に必要な最低カウント)。 |
| `DEFAULT_EXCLUDED_RESNAMES` | 既定の除外resname集合(結晶化添加物・イオン・修飾残基等、`frozenset[str]`)。 |

## `src/pymolrun`

PyMOLスクリプト(`.pml`)を専用conda/mamba環境で起動する。PyMOLは`rdkit>=2026.03.5`要件が共通環境
`pharmoforge`と競合するため専用環境(既定`pymol`)にインストールする前提(詳細は
[`alignview/README.md`](alignview/README.md#pymolの実行環境)参照)。当初`alignview`固有だったが
`pocketfinder`(`pf view-pocket`)でも同じ起動ロジックが必要になったためアトミックな技術要素として
切り出した。

### `pymolrun.launch`

| 関数 | 説明 |
| --- | --- |
| `run_pymol_script(script: str, pymol_env: str = "pymol") -> None` | `script`(PyMOLスクリプト本文)を一時`.pml`ファイルに書き出し、`mamba run -n <pymol_env> pymol <script>`でGUIモードのPyMOLを起動する。処理はPyMOLウィンドウを閉じるまでブロックし、実行後に一時ファイルを削除する。`mamba`/`conda`コマンドが見つからない場合は`RuntimeError`。 |

## `src/ligandprep`

SMILESからの3D配座生成・ドッキング用PDBQT変換([meeko](https://github.com/forlilab/meeko)、`pharmoforge`
環境にpipでインストール、`pyproject.toml`参照)。互変異性体・電荷状態の標準化は関知しない
(呼び出し側が[`molstd.standardize_smiles`](#srcmolstd)等で事前に行う前提)。使い方・環境構築は
[README.md](ligandprep/README.md)、実装経緯は[LIGANDPREP_PROMPT.md](ligandprep/LIGANDPREP_PROMPT.md)参照。

### `ligandprep.embed`

| 関数 | 説明 |
| --- | --- |
| `prepare_ligand_pdbqt(smiles: str, name: str, output_path: Path) -> Path` | SMILES1件から3D配座を1つ生成(RDKit ETKDGv3、固定シード)しMMFF94で最適化、meekoで電荷(Gasteiger)・原子タイプを割り当ててPDBQTとして`output_path`に書き出す。SMILESが不正、配座生成失敗、PDBQT変換失敗のいずれも`ValueError`。 |

## `src/docking`

受容体PDBQT準備(指定残基のフレキシブル化)・AutoDock Vinaの実行・結果パース・ポーズごとの受容体フル
コンフォメーション(PDB)+リガンドポーズ(SDF)の復元。vinaは`rdkit>=2026.03.5`とBoost.Pythonのビルドが
競合し`pharmoforge`環境に同居できないため、専用conda/mamba環境(既定`vina`、
[README.md](docking/README.md#実行環境)参照)に置く前提で`mamba run -n <env> vina ...`経由で起動する
(`pymolrun`と同じ回避パターン)。使い方は[README.md](docking/README.md)、実装経緯・処理内容の詳細は
[DOCKING_PROMPT.md](docking/DOCKING_PROMPT.md)参照。実データでの一連の使用例は
[cdk20_investigation.ipynb](notebooks/cdk20_investigation.ipynb)セクション7・7.1。

### `docking.receptor`

| 関数 | 説明 |
| --- | --- |
| `prepare_flexible_receptor(structure_path: Path, flexible_residues: list[tuple[str, int]], output_basename: Path) -> FlexReceptor` | 構造ファイル(PDB/CIF)から受容体PDBQTを準備する。`flexible_residues`((chain_id, resnum)のリスト)で指定した残基だけをmeeko(`Polymer.flexibilize_sidechain`)で可動側鎖として切り出し`<output_basename>_flex.pdbqt`に、残りは`<output_basename>_rigid.pdbqt`に書き出す(Vinaの`--flex`/`--receptor`にそれぞれ対応)。受容体全体のトポロジーを`<output_basename>.json`(`polymer.to_json()`)にも書き出す(`docking.export.export_docked_poses`がドッキング後の受容体フルコンフォメーション復元に使う)。水素付加・電荷割当(Gasteiger)はmeekoが内部で行う。指定残基が構造中に見つからない場合は`ValueError`。 |
| `FlexReceptor` | 準備結果(dataclass)。`rigid_pdbqt: Path`、`flex_pdbqt: Path \| None`(フレキシブル残基指定がなければ`None`)、`polymer_json: Path`、`n_flexible_residues: int`。 |

### `docking.vina`

| 関数 | 説明 |
| --- | --- |
| `calc_search_box(coords, padding: float = 4.0) -> tuple[tuple[float, float, float], tuple[float, float, float]]` | 座標配列(例: ポケット残基のCA座標)を包含するVina探索ボックスの中心・サイズ(Å)を計算する(meeko `gridbox.calc_box`のラッパー)。ポケット周辺残基のCAは既にポケットの縁まで広がっているため、paddingを大きくしすぎると探索空間が不必要に広がりドッキングが遅くなる。 |
| `run_vina(rigid_pdbqt: Path, ligand_pdbqt: Path, center, size, output_path: Path, flex_pdbqt: Path \| None = None, exhaustiveness: int = 8, num_modes: int = 9, cpu: int \| None = None, seed: int = 0, vina_env: str = "vina") -> VinaResult` | `mamba run -n <vina_env> vina --receptor ... --flex ... --ligand ... --out ...`でドッキングを実行し、`parse_vina_output`で出力ポーズのスコアを取得して返す。vinaが非0で終了、または`mamba`/`conda`コマンドが見つからない場合は`RuntimeError`。フレキシブル残基数が多いほど探索空間が急激に広がり実行時間が伸びる(実測: 3残基でリガンド1件約10秒、19残基では1件が3分超でも完了しない)。 |
| `parse_vina_output(output_path: Path) -> list[VinaPose]` | 出力PDBQTの`REMARK VINA RESULT:`行からポーズごとのスコアを抽出する(モード番号順)。`run_vina`が内部で使うほか、既に実行済みの出力(キャッシュ)を再読み込みする際にも呼べる。ポーズを1件も抽出できない場合は`ValueError`。 |
| `VinaResult` | ドッキング結果(dataclass)。`poses: list[VinaPose]`(モード番号順)、`output_path: Path`、`best_affinity: float`(プロパティ、`poses[0].affinity`)。 |
| `VinaPose` | 1ポーズ分のスコア(dataclass)。`mode: int`、`affinity: float`(kcal/mol)、`rmsd_lb: float`、`rmsd_ub: float`。 |

### `docking.export`

Vinaの`--out`出力自体にはリガンドポーズと可動側鎖(フレキシブル残基)の座標は含まれるが、受容体のリジッド
部分は含まれない(`--receptor`に渡した`_rigid.pdbqt`は不変のまま使い回されるため)。インタラクション解析
やMD初期構造として使える「ポーズごとの受容体フルコンフォメーション」を得るには、リジッド部分の構造
(`prepare_flexible_receptor`が書き出す`polymer_json`)とドッキング後の可動側鎖の座標(Vina出力)を
結合する必要があり、その結合をmeekoの`mk_export.py`(CLI)と同じ手順(`export_pdb_updated_flexres`)で行う。

| 関数 | 説明 |
| --- | --- |
| `export_docked_poses(polymer_json: Path, vina_output_pdbqt: Path, output_dir: Path, name: str, modes: list[int] \| None = None) -> list[ExportedPose]` | Vina出力(リガンド+フレキシブル受容体のPDBQT)から、モードごとに受容体のフルコンフォメーション(リジッド部分+可動側鎖、標準PDB形式・水素付き)を`<output_dir>/<name>_mode<N>_receptor.pdb`に、リガンドポーズ(結合次数を復元したSDF)を`<output_dir>/<name>_mode<N>_ligand.sdf`に書き出す。`modes`省略時は全モードを書き出す(1始まり、Vinaのスコア順)。SDF変換に失敗した場合は`ValueError`。 |
| `ExportedPose` | 1ポーズ分の書き出し結果(dataclass)。`mode: int`、`receptor_pdb: Path`、`ligand_sdf: Path`。 |

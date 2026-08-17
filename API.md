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
| `resolve_uniprot_accession(identifier: str) -> str` | UniProt entry name / accessionのいずれを与えてもUniProt accessionを返す(`fetcher`の`structure=` `--af`指定時、ダウンロード用accessionと出力ファイル名(`<accession>_AF.<fmt>`)の両方の解決に使用)。 |
| `pdb_id_to_uniprot_accessions(pdb_id: str) -> list[str]` | PDB ID配下の全ポリマーエンティティに紐づくUniProt accessionを重複なく出現順で返す(RCSB Data API `data.rcsb.org`を使用)。複合体で複数の蛋白質が含まれる場合は複数返る。`fetcher`の`structure=` `--type=fasta`(`--af`なし)時、PDB IDからUniProt accessionを解決するのに使用。 |

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
| `FitResult` | `fit_by_residue_number`の戻り値(dataclass)。`matrix`(4x4 numpy配列、行優先。`v' = matrix @ [x, y, z, 1]`)、`rmsd: float`、`n_residues: int`、`mobile_chain: str`、`target_chain: str`。 |

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

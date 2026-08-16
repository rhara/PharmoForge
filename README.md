# PharmoForge

創薬支援ツール群のプロジェクト。

## 想定機能一覧

現時点での想定であり、まだ全機能を網羅しているわけではない。これは大まかな機能分類であり、実際の`src`以下のパッケージ構成はこれより技術的な観点で細かく分割される([パッケージ構成](#パッケージ構成)、[PharmoForgeのCLAUDE.md](CLAUDE.md)参照)。

| 機能 | パッケージ名(目安) |
| --- | --- |
| データ(物質、実験データ、構造、科学的情報)の収集 | `fetcher` |
| 活性データの探索的解析(scaffold比較等のSAR分析) | `scaffoldanalyzer` |
| 活性予測(機械学習) | `activitypredictor` |
| 物性予測(機械学習) | `propertypredictor` |
| 蛋白構造の準備 | `proteinprep` |
| 蛋白のシーケンス・構造解析 | `proteinanalyzer` |
| リガンドの構造準備(コンフォメーション、トートマー、標準化、プロトネーション等) | `ligandprep` |
| ドッキング(リジッド、フレキシブル、コバレント等) | `docking` |
| フォルディング予測 | `folding` |
| 分子動力学シミュレーション(metadynamics, steered, replica exchange等種々) | `dynamics` |
| 分子動力学結果の解析 | `dynamicsanalizer` |
| 分子動力学の応用(free energy関連、リガンド安定性評価、FEP、MM/XBSA等) | `freeenergy` |

## パッケージ構成

上記の機能一覧は大まかな分類であり、アトミックで独立してテスト・差し替えがしやすい単位は、現に複数の機能から使われているかに関わらず、保守性のため専用の小さなパッケージに切り出す方針。

現時点で実装済みのパッケージ:

| パッケージ | 役割 |
| --- | --- |
| `src/core` | 汎用ユーティリティ(verboseなログ出力等) |
| `src/idmap` | 蛋白識別子マッピング(UniProt entry name / accession / ChEMBL target id相互変換)。当初`fetcher`固有だったがアトミックな技術要素として独立させた |
| `src/chembl` | ChEMBL REST APIからの活性データ取得。当初`fetcher`固有だったがアトミックな技術要素として独立させた |
| `src/molstd` | 化合物構造の標準化(RDKit標準の`rdMolStandardize`のみを用いた自前実装、外部パッケージ`chembl_structure_pipeline`には非依存)。当初`fetcher`固有だったがアトミックな技術要素として独立させた |
| `src/rcsb` | RCSB PDBからの構造ファイルダウンロード。当初`fetcher`固有だったがアトミックな技術要素として独立させた |
| `src/afdb` | AlphaFold DBからの構造ファイルダウンロード。当初`fetcher`固有だったがアトミックな技術要素として独立させた |
| `src/uniprot` | UniProtエントリの取得・創薬向け情報抽出。当初`proteinanalyzer`固有だったがアトミックな技術要素として独立させた |
| `src/molscaffold` | 化合物のBemis-Murckoスキャフォールド計算。当初`scaffoldanalyzer`固有だったがアトミックな技術要素として独立させた |
| `src/actbin` | 活性値の分位点ビニング(high/mid/low分類)。当初`scaffoldanalyzer`固有だったがアトミックな技術要素として独立させた |
| `src/structfit` | 残基番号の対応のみに基づく剛体重ね合わせ計算(ProDy)。`alignview`(`--method number`)向けにアトミックな技術要素として切り出した |
| `src/structio` | PDB/CIF構造ファイルの読み書き(拡張子で自動判別、ProDy)。`structfit`/`proteinextract`の共通処理としてアトミックな技術要素に切り出した |
| `src/fetcher` | ChEMBL活性データの標準化・集約・TSV書き出し、RCSB PDB/AlphaFold DB構造データ取得のCLI([詳細](fetcher/README.md)) |
| `src/scaffoldanalyzer` | Bemis-Murckoスキャフォールド単位での活性高低の分布比較([詳細](scaffoldanalyzer/README.md)) |
| `src/proteinprep` | PDBFixerによる蛋白構造の欠損原子補完・プロトン化([詳細](proteinprep/README.md)) |
| `src/proteinanalyzer` | UniProtからの創薬向け蛋白情報取得CLI([詳細](proteinanalyzer/README.md)) |
| `src/proteinextract` | 構造ファイルから指定チェーンの抽出・水分子の除去([詳細](proteinextract/README.md)) |
| `src/alignview` | 複数のPDB/CIF構造をPyMOLで開き先頭構造にアラインして表示([詳細](alignview/README.md)) |

`pf`コマンド自体は`src/core/cli.py`のclickグループが起点となり、各機能パッケージがサブコマンドを登録する。
サブコマンド一覧は[`pf`コマンド一覧](#pfコマンド一覧)を参照。

アトミックなパッケージ(`core`/`idmap`/`chembl`/`molstd`/`rcsb`/`afdb`/`uniprot`/`molscaffold`/`actbin`/`structfit`/`structio`等)の関数一覧は[API.md](API.md)を参照。

## `pf`コマンド一覧

各コマンドの詳細な引数・出力仕様は、それぞれのリンク先の機能READMEを参照。

| コマンド | 機能 | 例 |
| --- | --- | --- |
| [`pf fetch`](fetcher/README.md) | ChEMBL活性データ、RCSB PDB/AlphaFold DB構造データの取得 | `pf fetch activity=CDK4_HUMAN --outdir data` |
| [`pf analyze-scaffolds`](scaffoldanalyzer/README.md) | Bemis-Murckoスキャフォールド単位での活性高低分布の比較 | `pf analyze-scaffolds data/cdk4_human_activities.tsv --output-dir data/cdk4_scaffold_analysis` |
| [`pf prep-protein`](proteinprep/README.md) | PDBFixerによる蛋白構造の欠損原子補完・プロトン化 | `pf prep-protein data/9csk.cif --output data/9csk_dock.pdb --mode dock` |
| [`pf protein-info`](proteinanalyzer/README.md) | UniProtからの創薬向け蛋白情報取得(JSON出力) | `pf protein-info EGFR_HUMAN --output data/egfr_info.json` |
| [`pf protein-extract`](proteinextract/README.md) | 構造ファイルから指定チェーンの抽出・水分子の除去 | `pf protein-extract data/cdk2/2CCH.cif --chains=A,G,H,I,J --remove-water --output data/cdk2/2CCH_main.cif` |
| [`pf align-view`](alignview/README.md) | 複数のPDB/CIF構造をPyMOLで開き先頭構造にアラインして表示 | `pf align-view data/tyk2/TYK2_HUMAN_af.cif data/tyk2/6NZP.cif data/tyk2/4OLI.cif data/tyk2/5C03.cif` |

全コマンド共通で `pf <コマンド> --help` により引数の詳細を確認できる。

## 依存パッケージ

- パッケージは可能な限りconda-forgeからインストールする。
- RDKitは中核をなすライブラリの一つであり、`rdkit >= 2026.03.5` でpinする。
- これを満たせないパッケージのインストールが必要な場合は、そのパッケージ専用の別環境を作る(該当数は多くないはず)。

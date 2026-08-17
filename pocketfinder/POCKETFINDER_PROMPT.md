# pocketfinder 実装記録

このドキュメントは `pocketfinder` 機能を再現するための仕様記録。

## 目的

蛋白のポケット(リガンド結合部位候補)を検出し、各ポケットに面する残基をリストアップする。
PharmoForgeにおける「蛋白のシーケンス・構造解析」機能の一つ(ドッキング前のポケット同定、
変異設計や結合部位比較の基礎データとして利用)。

## 手法の選定

ポケット検出には[fpocket](https://github.com/Discngine/fpocket)(Voronoi分割ベースの幾何学的手法)を
採用した。理由:

- conda-forgeから配布されており(`fpocket`パッケージ)、GPU非搭載環境でもそのまま動作する
  (`CLAUDE.md`のGPU非搭載環境でのCPUフォールバック要件を満たす。本手法はそもそもGPUを使わない)。
- P2Rank等の機械学習ベース手法と比べ依存が軽く、追加のモデル取得が不要。
- 出力にポケットスコア・druggabilityスコア・ポケットに面する原子(≒残基)一覧が標準で含まれ、
  「ポケット周辺残基のリストアップ」という要件にそのまま対応する。

## fpocketのCLI仕様(実データで確認)

```bash
fpocket -f <入力PDB/CIFファイル> -w pdb
```

- 入力ファイルと同じディレクトリに`<stem>_out/`(`stem`は入力ファイル名から拡張子を除いた部分)を
  生成する。`-o`等の出力先指定オプションは存在しない。
- `<stem>_out/<stem>_info.txt`: 検出された各ポケット(`Pocket N :`ブロック)のスコア・druggability
  スコア・alphaスフィア数・SASA・体積等の記述子。フォーマットに`Key : \tValue`と`Key:\tValue`が
  混在するため、パースは最初の`:`で分割し両辺をstripする方式で統一的に扱う
  (`src/pocket/fpocket.py`の`_parse_info_file()`)。
- `<stem>_out/pockets/pocket<N>_atm.<拡張子>`: ポケットNに面する原子(ポケットの
  Voronoi頂点(alphaスフィア)に接触する蛋白原子)。既定では入力と同じ形式(pdb/cif)で出力されるが、
  fpocketが生成するmmCIFは`pdbx_PDB_model_num`列を欠くためProDy(`structio.parse_structure`)で
  パースできないことを実データで確認した。このため`-w pdb`を常に指定し、ポケット周辺残基抽出用の
  出力は常にPDB形式に固定している(`<stem>_out.<形式>`本体やPyMOL/VMDスクリプトは入力と同じ形式で
  出力される)。
- ポケットは`Pocket 1`がスコア最高、以降スコア降順で並ぶことを実データで確認済みだが、念のため
  `run_fpocket()`側でも`score`降順に再ソートしている。

## 実装ファイル

- `src/pocket/fpocket.py` — fpocketの実行・出力パース(アトミックな技術要素、`API.md`参照)。
  - `run_fpocket(structure_path, work_dir) -> list[Pocket]`: `structure_path`を`work_dir`にコピーし
    `fpocket -f <コピー> -w pdb`を実行、`<stem>_info.txt`と`pockets/pocket<N>_atm.pdb`をパースして
    `Pocket`(`pocket_id`/`score`/`druggability_score`/`n_alpha_spheres`/`volume`/`total_sasa`/
    `polar_sasa`/`apolar_sasa`/`hydrophobicity_score`/`residues: list[PocketResidue]`)のリストを
    スコア降順で返す。`PocketResidue`は`chain_id`/`resnum`/`resname`(重複除去、chain/resnum昇順)。
    fpocketが非0で終了した場合は`RuntimeError`。
- `src/pocketfinder/report.py` — `list[Pocket]`のJSON/TSV整形・書き出し・読み込み(pocketfinder固有)。
  - `format_pockets_json()`/`write_pockets_json()`: `pf find-pocket`の出力用。
  - `read_pockets_json()`: `pf view-pocket`/`pf show-pocket`が`pockets.json`を`list[Pocket]`に
    復元するための読み込み。
  - `pockets_to_dataframe()`: `list[Pocket]`を1行1ポケットの`pandas.DataFrame`に変換する
    (`residues`列は`<chain>:<resnum>`をカンマ結合した1セルにまとめる)。
  - `format_pockets_table()`/`write_pockets_table()`: 上記DataFrameを`to_csv(sep="\t")`でTSVに
    整形・書き出し(`pf show-pocket`の出力用、[`scaffoldanalyzer.report.write_summary_tsv()`]
    (../scaffoldanalyzer/README.md)と同じ方式)。
- `src/pocketfinder/view.py` — ポケット可視化用PyMOLスクリプトの組み立てと起動(pocketfinder固有)。
  - `build_pocket_view_script(structure_path, pockets, top_n=None) -> str`: 構造を読み込み、
    ポケットごとに周辺残基を`select`し、色分け(スティック表示+半透明サーフェス)するPyMOL
    スクリプト本文を組み立てる。残基選択はチェーンごとにグルーピングした
    `(chain <ID> and resi <番号+番号+...>)`のOR結合。オブジェクト名を構造ファイル名(`protein`等の
    予約語)ではなくファイル名(拡張子なし)にする(`protein`は`polymer.protein`セレクタと衝突し警告が
    出ることを実データで確認)。
  - `launch_pocket_view(structure_path, pockets, pymol_env="pymol", top_n=None) -> None`:
    スクリプトを組み立て[`pymolrun.run_pymol_script()`](../src/pymolrun)で起動する。
- `src/pocketfinder/cli.py` — `pf find-pocket`/`pf view-pocket`/`pf show-pocket`サブコマンド。

`pf show-pocket`の行単位はポケット単位(1行1ポケット)とした(残基単位ではなく)。理由:
ポケット同士のスコア・druggability等の比較が主目的であり、ポケット単位の方が一覧性が高い
(ユーザーとの選択式確認で決定)。残基一覧は`residues`列に`<chain>:<resnum>`のカンマ結合で
まとめて含める(TSVはタブ区切りのためカンマは値の一部として問題なく扱える)。

## `src/pymolrun`への切り出し

PyMOLスクリプトを専用conda/mamba環境で起動する処理(`mamba run -n <env> pymol <script>`、一時
`.pml`ファイルの書き出し・実行後の削除)は、元々`alignview/view.py`の`launch_alignment_view()`に
実装されていた。`view-pocket`でも同じ起動ロジックが必要になったため、`run_pymol_script(script,
pymol_env="pymol")`として`src/pymolrun/launch.py`にアトミックな技術要素として切り出し、
`alignview.view.launch_alignment_view()`もこれを使うようリファクタした(挙動は変更なし)。

## CLI仕様

```
pf find-pocket <構造ファイル> --outdir <出力ディレクトリ> [--top N]
pf view-pocket <構造ファイル> --pockets <pockets.jsonのパス> [--top N] [--pymol-env <env名>]
pf show-pocket <pockets.jsonのパス> [--output <出力TSVファイル>]
```

`find-pocket`の処理の流れ(`src/pocketfinder/cli.py`):

1. `pocket.run_fpocket(structure_path, output_dir)`でポケット検出・残基抽出を実行
   (fpocketの生出力一式も`output_dir`に残る)。
2. `--top`指定時は上位N件に絞り込み。
3. `src/pocketfinder/report.write_pockets_json()`で`output_dir/pockets.json`に書き出し。

`view-pocket`の処理の流れ:

1. `--pockets`で指定された`pockets.json`を`report.read_pockets_json()`で`list[Pocket]`に復元。
2. `view.launch_pocket_view(structure_path, pockets, pymol_env=..., top_n=...)`でPyMOLを起動
   (`--top`はここで渡され、スクリプト組み立て時に絞り込まれる)。

`show-pocket`の処理の流れ:

1. 位置引数で指定した`pockets.json`を`report.read_pockets_json()`で`list[Pocket]`に復元
   (`find-pocket`/`view-pocket`と異なり、構造ファイル自体は不要なため引数に取らない)。
2. `--output`省略時は`report.format_pockets_table()`で整形したTSVを`click.echo(..., nl=False)`で
   標準出力へ(`to_csv()`が末尾改行を含むため、`click.echo`側での改行付与は抑制)。指定時は
   `report.write_pockets_table()`でファイルに書き出し。

## PyMOLの実行環境

`view-pocket`は`alignview`と同様、専用conda/mamba環境(既定`pymol`)にインストールした
`pymol-open-source`を`mamba run -n <env> pymol <script>`経由で起動する前提
([alignview/README.md](../alignview/README.md#pymolの実行環境)参照)。

```bash
mamba create -n pymol -c conda-forge pymol-open-source
```

## テスト

```bash
pytest tests/pocketfinder tests/pocket tests/pymolrun
```

`tests/pocket/test_fpocket.py`は`_parse_info_file()`/`_extract_residues()`を実データ相当の合成
テキスト・PDBで実処理検証し、`run_fpocket()`は`subprocess.run`を`unittest.mock`でモック
(モック内で`<stem>_out/`の出力ファイルを合成生成してパース処理を検証)している。実際の`fpocket`
バイナリ呼び出しはテストでは行わない。`tests/pocketfinder/test_view.py`は`build_pocket_view_script()`
を実処理検証し、`launch_pocket_view()`は`pymolrun.run_pymol_script`をモックする。実際のPyMOL起動は
テストでは行わないが、生成したスクリプトを実際に`pymol -cq`(GUIなしのheadlessモード)で実行し、
エラーなく完了することを実データ(CDK2 AlphaFold構造 + 検出済みポケット3件)で確認済み。
`tests/pocketfinder/test_report.py`の`format_pockets_table()`/`write_pockets_table()`のテストは、
出力したTSVを`pandas.read_csv(sep="\t")`で読み戻して行数・列値を検証する。

## 動作例(サンプルデータ)

CDK2(ヒト、AlphaFold予測構造、UniProt: P24941)を題材にした検出例(実際に動作確認済み。
18ポケットを検出、最上位ポケット(pocket 1、score 0.188)はGly-richループ(14-15番残基)や
ATP結合部位近傍の残基(126-129番、149番、154-188番付近)を含む):

```bash
pf find-pocket data/cdk2/P24941_AF.cif --outdir data/cdk2/P24941_AF_pockets --top 3
pf view-pocket data/cdk2/P24941_AF.cif --pockets data/cdk2/P24941_AF_pockets/pockets.json --top 3
pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json | column -t
```

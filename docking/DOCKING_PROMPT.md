# docking 実装記録

このドキュメントは `docking` 機能を再現するための仕様記録。関数シグネチャ・挙動の正は
[API.md](../API.md#srcdocking)を参照。

## 目的

指定残基をフレキシブル(可動側鎖)にした状態でAutoDock Vinaドッキングを行い、スコアだけでなく
蛋白コンフォメーションとリガンドポーズをセットで(後のインタラクション解析・MD初期構造としての
再利用に耐える形で)保存する。PharmoForgeにおける「ドッキング」機能を担う
(現時点はフレキシブルドッキングのみ、リジッド・コバレントドッキングは未実装)。

現時点は最小限のスコープで実装している。以下は今後の拡張候補で未実装:

- `pf`コマンド化(現在は`pharmoforge`ライブラリ関数として直接呼び出す運用、[README.md](../README.md)参照)
- コバレントドッキング
- 複数リガンド・複数レプリカのバッチ実行の並列化(現状は逐次ループ、ノートブック側で実装)。
  Uni-Dockの`paired_batch`モード(1ポケット×多数リガンドのGPUバッチ処理向け)は理論上この用途に
  適合しうるが、フレキシブル残基指定(`--flex`)との併用可否が未検証(下記「エンジン選定」参照)

## 経緯: vinaを専用conda環境に分離した理由

AutoDock Vinaのconda-forgeビルド(`vina`パッケージ)は特定バージョンのlibboost-pythonを要求し、これが
`pharmoforge`環境のrdkit(`rdkit>=2026.03.5`、これも同じlibboost-pythonに依存)と競合する。実際に
`pharmoforge`環境へ`vina`をインストールしたところ、依存解決の結果rdkitが2026.03.5→2025.09.6へ強制
ダウングレードされることを確認した(pyproject.tomlの`rdkit>=2026.03.5`要件に抵触)。

CLAUDE.mdの「これを満たせないパッケージのインストールが必要な場合は、そのパッケージ専用の別環境を作る」
という既存方針に従い、vinaは専用env(`vina`、[README.md](../README.md#実行環境)参照)にインストールし、
`pymolrun.run_pymol_script`と同じパターン(`mamba run -n <env> vina ...`のsubprocess起動)で呼び出す
こととした。

## 経緯: meeko/gemmiをpipインストールした理由

`meeko`(受容体・リガンドのPDBQT準備、ポーズ復元に使用)とその依存`gemmi`のconda-forgeビルドは
python<3.14までにしか対応しておらず、`pharmoforge`環境(python 3.14)にconda-forgeから直接
インストールできない。両者ともPyPIのwheelはpure Python(`meeko`)/cp314向けmanylinux wheel(`gemmi`)が
存在し、python 3.14での動作を確認済みのため、CLAUDE.mdの「これを満たせないパッケージ...のみpip許可」
に従いpipでインストールする(`pyproject.toml`にコメント付きで記載)。

## エンジン選定(smina/gninaとの比較)

フレキシブル残基数が多いとVinaのドッキングが遅くなる問題(下記「実行時間について」参照)に対し、
smina・gnina(GPU)を同一条件で実測比較した。CDK20 AlphaFold構造・フレキシブル残基15個・
`exhaustiveness=8`・`cpu=4`・同一リガンドでの結果:

| エンジン | 所要時間 | 備考 |
| --- | --- | --- |
| vina | 707秒 | CPU |
| smina | 805秒(vinaより14%遅い) | CPU |
| gnina(`--scoring vina --cnn_scoring none`) | 736秒 | GPU使用(最小化ステップでGPU使用率41%を確認)だが有意な高速化なし |

smina・gninaはいずれもvinaのフォークで、探索アルゴリズム(Monte Carlo + 局所最適化)そのものは共通。
遅さの原因は「フレキシブル残基の回転可能な結合数に応じて組合せ的に増大する探索空間」であり、
GPU化されている部分(最小化の数値計算等)を高速化しても支配的なボトルネックにはならないため、
gninaでGPUが実際に稼働していても全体時間はvinaとほぼ変わらなかった。

Uni-Dock(GPU、公称vina比1000倍以上)も調査したが、公式ドキュメントで「単一リガンドのフレキシブル
ドッキングはオーバーヘッドの割合が大きく、かなり遅くなる」と明記されており、高速化は「1ポケットに
対し大量(1000規模)のリガンド」をバッチ処理する用途向け。本ノートブックの実際のワークフロー
(1ポケット×化合物40件)はこの形に近く理論上は有望だが、バッチモード(`paired_batch`)と
`--flex`(フレキシブル残基指定)を併用できるかどうかが公式ドキュメントに記載がなく未検証のため、
今回は導入を見送った(実際に試すには専用env構築+実データでの動作確認が必要、上記「今後の拡張候補」)。

以上の実測・調査に基づき、現時点ではエンジンをvinaから変更する明確なメリットがないと判断し、
vinaを継続採用している。フレキシブル残基数・exhaustivenessの調整で対応する方針(下記「実行時間に
ついて」参照)。

## CLI仕様

なし。`pf`コマンド化は現在保留中、Pythonから直接呼び出す(使い方は[README.md](../README.md#使い方)参照)。

### 処理内容

#### 1. 受容体準備(`src/docking/receptor.py` の `prepare_flexible_receptor()`)

1. 構造ファイル(PDB/CIF)をProDy(`prody.parsePDB`/`prody.parseMMCIF`)で読み込む。
2. `ResidueChemTemplates.create_from_defaults()` + `MoleculePreparation()` で
   `Polymer.from_prody(...)` を構築する(meekoの受容体トポロジー表現)。
3. 指定残基((chain_id, resnum)のリスト)ごとに `polymer.flexibilize_sidechain(res_id, mk_prep)` を
   呼び、可動側鎖として切り出す。指定残基が受容体に存在しない場合は`ValueError`。
   側鎖に回転可能な結合を持たない残基(Gly/Ala)を指定した場合、meeko側が「no movable atoms」の
   警告を出しつつ静かにスキップする(呼び出し側でGly/Alaを事前に除外することを推奨、
   [cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション7参照)。
4. `PDBQTWriterLegacy.write_from_polymer(polymer)` でリジッド部分(`<basename>_rigid.pdbqt`)と
   フレキシブル残基ごとのPDBQT(結合して`<basename>_flex.pdbqt`)を書き出す
   (フレキシブル指定が0件、または全て可動原子なしだった場合は`_flex.pdbqt`を書き出さず`None`)。
5. `polymer.to_json()` を`<basename>.json`にも書き出す。Vinaの`--receptor`/`--flex`入力自体には
   リジッド部分の座標が(ドッキング中不変のため)含まれず、ドッキング後に受容体フルコンフォメーションを
   復元する際(手順3「ポーズの復元」)にこのJSONが必要になる。

#### 2. 探索ボックス(`src/docking/vina.py` の `calc_search_box()`)

meeko `gridbox.calc_box(coords, padding)` のラッパー。座標配列(min/max)を包含する中心・サイズを
計算する。既定`padding=4.0`(Å)。ポケット周辺残基のCA座標をそのまま渡す運用を想定しており、それらは
既にポケットの縁まで広がっているため、paddingを大きくしすぎない(実測で8.0だと1辺30Å超になり
探索が不必要に遅くなることを確認、詳細は「実行時間について」参照)。

#### 3. ドッキング実行(`src/docking/vina.py` の `run_vina()`)

1. `shutil.which("mamba")`(なければ`conda`)を探す。見つからなければ`RuntimeError`。
2. `mamba run -n <vina_env> vina --receptor <rigid_pdbqt> [--flex <flex_pdbqt>] --ligand <ligand_pdbqt>
   --center_x/y/z ... --size_x/y/z ... --out <output_path> --exhaustiveness <N> --num_modes <N>
   --seed <N> [--cpu <N>]` をsubprocess実行する。非0終了で`RuntimeError`。
3. `parse_vina_output(output_path)` で結果をパースして返す。

#### 4. 結果パース(`src/docking/vina.py` の `parse_vina_output()`)

出力PDBQTの`REMARK VINA RESULT:\s+(affinity)\s+(rmsd_lb)\s+(rmsd_ub)`行を正規表現で抽出し、
モード番号順(Vinaのスコア順)の`VinaPose`リストを返す。1件も見つからない場合は`ValueError`。
`run_vina`が内部で使うほか、既存の出力ファイルをドッキングし直さずに再集計する際にも呼べる
(ノートブック側のキャッシュ機構で使用)。

#### 5. ポーズの復元(`src/docking/export.py` の `export_docked_poses()`)

meekoの`mk_export.py`(CLI)と同じ手順を踏む:

1. `Polymer.from_json(...)` で手順1が書き出した`.json`を読み込む(リジッド部分を含む受容体全体の
   トポロジー)。
2. `PDBQTMolecule.from_file(vina_output_pdbqt, skip_typing=True)` でVina出力(全モード)を読み込む。
3. `RDKitMolCreate.write_sd_string(pdbqt_mol, keep_flexres=False)` でリガンド側を、結合次数を
   復元したSDF文字列(複数モード分、`$$$$`区切り)に変換する。空文字列(=全モード変換失敗)は
   `ValueError`。
4. モードごとに、`pdbqt_mol`を1ポーズ分に絞った複製(`_positions`/`_pose_data["n_poses"]`/
   `_current_pose`を書き換え)を作り、`export_pdb_updated_flexres(polymer, single_pose)` で
   受容体側を、そのモードの可動側鎖座標で更新した標準PDB文字列(リジッド部分+水素付き)として
   復元する(`polymer`・`pdbqt_mol`とも呼び出しごとに複製し、副作用を防ぐ)。
5. `<output_dir>/<name>_mode<N>_receptor.pdb` / `<output_dir>/<name>_mode<N>_ligand.sdf` として
   書き出す。`modes`引数(省略時は全モード)で書き出すモードを絞れる。

実データ(CDK20 AlphaFold構造、PHE81をフレキシブル指定)で、モード間で側鎖原子(CD1/CE1/CZ等)の
座標が実際に異なる(=ドッキングごとに正しく更新されている)ことを確認済み。

## 実装ファイル

- `src/docking/receptor.py` — フレキシブル受容体PDBQT準備
- `src/docking/vina.py` — 探索ボックス計算・Vina実行・結果パース
- `src/docking/export.py` — ポーズごとの受容体フルコンフォメーション・リガンドポーズの復元

リガンド側のPDBQT準備は本機能の対象外。[`ligandprep`](../ligandprep/LIGANDPREP_PROMPT.md)機能を使う。

## 依存パッケージ

- `vina`(AutoDock Vina、専用env`vina`、conda-forge)。理由は上記「経緯」参照。
- `meeko`/`gemmi`(`pharmoforge`環境にpip)。理由は上記「経緯」参照。

## テスト

ネットワークアクセスを伴わないため、実データ(最小限の合成PDBフラグメント・小分子SMILES)を使った
実処理での検証を基本とする。Vinaのsubprocess呼び出し自体はモックする(`docking.vina`のテスト、
[`pymolrun`](../src/pymolrun)のテストパターンを踏襲)。

```bash
pytest tests/docking tests/ligandprep
```

## 動作例(サンプルデータ)

[cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション7・7.1:
fpocket + 保存モチーフ + 相同蛋白の共結晶化リガンド接触残基から特定したATP結合部位周辺残基
(Gly/Ala除く)をフレキシブルにし、ChEMBLから収集した「他のCDKで高活性を示した化合物」を
CDK20のAlphaFold予測構造にドッキングし、各化合物の最良ポーズについて受容体コンフォメーション(PDB)
とリガンドポーズ(SDF)を書き出す一連の流れ。

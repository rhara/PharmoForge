# docking

指定残基をフレキシブル(可動側鎖)にした[AutoDock Vina](https://vina.scripps.edu/)ドッキングを行う機能。
受容体PDBQT準備・Vina実行・結果パースに加え、ドッキング後の受容体フルコンフォメーション(PDB)とリガンド
ポーズ(SDF)をポーズごとに復元する(インタラクション解析・MD初期構造としての利用を想定)。

`pf`コマンド化は現在保留中(トップ[README.md](../README.md)参照)。現時点ではPythonから直接呼び出す。
リガンド側の入力は[`ligandprep`](../ligandprep/README.md)で用意する。
関数シグネチャの正は[API.md](../API.md#srcdocking)を参照(ここでは典型的な使い方の流れのみ示す)。
実データでの使用例は[cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション7・7.1。

## 実行環境

vinaは`rdkit>=2026.03.5`とBoost.Pythonのビルドが競合し`pharmoforge`環境に同居できないため、専用の
conda/mamba環境が必要。

```bash
mamba create -n vina -c conda-forge python=3.14 vina
```

(`docking.run_vina`の`vina_env`引数、既定`vina`)

## 使い方

```python
from ligandprep import prepare_ligand_pdbqt
from docking import prepare_flexible_receptor, calc_search_box, run_vina, export_docked_poses

# 1. 受容体準備: 指定残基((chain_id, resnum)のリスト)をフレキシブルにする
flex_receptor = prepare_flexible_receptor(
    "data/cdk20/Q8IZL9_repaired.pdb",
    [("A", 10), ("A", 33), ("A", 81), ("A", 82)],
    "data/cdk20/docking/cdk20",
)

# 2. 探索ボックス: ポケット周辺残基のCA座標を包含する範囲を計算
box_center, box_size = calc_search_box(pocket_ca_coords)

# 3. リガンド準備(ligandprep)
ligand_pdbqt = prepare_ligand_pdbqt("CC(=O)Nc1ccc(O)cc1", "acetaminophen", "data/cdk20/ligands/acetaminophen.pdbqt")

# 4. ドッキング実行
result = run_vina(
    rigid_pdbqt=flex_receptor.rigid_pdbqt,
    ligand_pdbqt=ligand_pdbqt,
    center=box_center,
    size=box_size,
    output_path="data/cdk20/poses/acetaminophen_docked.pdbqt",
    flex_pdbqt=flex_receptor.flex_pdbqt,
)
print(result.best_affinity)  # kcal/mol

# 5. ポーズごとの受容体コンフォメーション+リガンドポーズを復元(最良ポーズのみ、全モードはmodes=None)
exported = export_docked_poses(
    polymer_json=flex_receptor.polymer_json,
    vina_output_pdbqt=result.output_path,
    output_dir="data/cdk20/exported_poses",
    name="acetaminophen",
    modes=[1],
)
print(exported[0].receptor_pdb, exported[0].ligand_sdf)
```

実例は[cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション7を参照。

### 処理内容

1. **受容体準備**(`prepare_flexible_receptor`): 構造ファイル(PDB/CIF)から、指定残基だけを可動側鎖と
   して切り出したPDBQT(`_rigid.pdbqt`/`_flex.pdbqt`)を書き出す。受容体全体のトポロジー(`.json`)も
   書き出す(手順5で使う)。
2. **探索ボックス**(`calc_search_box`): 座標配列(通常はポケット周辺残基のCA)を包含するボックスの
   中心・サイズを計算する(meeko `gridbox`のラッパー)。paddingを大きくしすぎると探索空間が不必要に
   広がりドッキングが遅くなるため、既定は控えめな4Å。
3. **ドッキング実行**(`run_vina`): 専用env(既定`vina`)で`vina`実行ファイルをsubprocess起動する。
4. **結果パース**(`parse_vina_output`): 出力PDBQTの`REMARK VINA RESULT:`行からポーズごとのスコアを
   抽出する。キャッシュされた既存の出力を再読み込みする際にも使える。
5. **ポーズの復元**(`export_docked_poses`): Vina出力自体には可動側鎖とリガンドの座標しか含まれず、
   受容体のリジッド部分(不変)は含まれない。手順1の`.json`(受容体全体のトポロジー)と組み合わせて、
   ポーズごとの受容体フルコンフォメーション(標準PDB形式・水素付き)とリガンドポーズ(結合次数を
   復元したSDF)を書き出す。

## 実行時間について

フレキシブル残基数が探索空間・実行時間に大きく効く。CDK20 AlphaFold構造での実測: hinge付近の3残基のみ
フレキシブルにした場合はリガンド1件・`exhaustiveness=8`で約10秒。ATP結合部位周辺の19残基(Gly/Ala除く)
全てをフレキシブルにすると1件が3分経っても完了しない。フレキシブル残基数・`exhaustiveness`・化合物数は
用途に応じて調整すること。大規模バッチではキャッシュ(既存ファイルはスキップ)を前提にした分割実行・
中断再開を推奨する(ノートブックの実装を参照)。

## 実装方針

- ドッキングエンジンは[AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)(専用env、上記参照)。
- 受容体PDBQT準備・ポーズ復元は[meeko](https://github.com/forlilab/meeko)を使う(`ligandprep`と同じ理由で
  `pharmoforge`環境にpipインストール)。
- GPU対応(AutoDock-GPU/gnina等)は未実装。現状はVinaのマルチスレッド(`--cpu`)のみ。

## テスト

```bash
pytest tests/docking
```

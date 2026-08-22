# ligandprep

SMILESからの3D配座生成・ドッキング用PDBQT変換を行う機能。互変異性体・電荷状態の標準化は関知しない
(呼び出し側が[`molstd`](../src/molstd)等で事前に行う前提)。`docking`のリガンド側入力を用意するために使う。

`pf`コマンド化は現在保留中(トップ[README.md](../README.md)参照)。現時点ではPythonから直接呼び出す。
関数シグネチャの正は[API.md](../API.md#srcligandprep)を参照(ここでは典型的な使い方のみ示す)。
実データでの使用例は[cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション7。

## 使い方

```python
from ligandprep import prepare_ligand_pdbqt

prepare_ligand_pdbqt("CC(=O)Nc1ccc(O)cc1", "acetaminophen", "data/ligands/acetaminophen.pdbqt")
```

### 処理内容

1. RDKitでSMILESをパースし、明示的水素を付加する。
2. ETKDGv3(固定シード)で3D配座を1つ生成し、MMFF94で最適化する。
3. meeko(`MoleculePreparation`)でGasteiger電荷・原子タイプを割り当て、PDBQTとして書き出す。

配座生成・PDBQT変換に失敗した場合は`ValueError`。

## 実装方針

- 3D配座生成・構造最適化はRDKitのみで行う(外部ツール非依存)。
- PDBQT変換は[meeko](https://github.com/forlilab/meeko)を使う。conda-forgeのビルドがpython<3.14までにしか
  対応していないため、`pharmoforge`環境にpipでインストールする(詳細は[README.md](../README.md#依存パッケージ)参照)。

## テスト

```bash
pytest tests/ligandprep
```

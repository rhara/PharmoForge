# ligandprep 実装記録

このドキュメントは `ligandprep` 機能を再現するための仕様記録。関数シグネチャ・挙動の正は
[API.md](../API.md#srcligandprep)を参照。

## 目的

SMILES(標準化済み想定)1件から、ドッキングにそのまま使えるPDBQTファイルを生成する。
PharmoForgeにおける「リガンドの構造準備」機能を担う。`docking`機能のリガンド側入力を用意するために作った。

互変異性体・電荷状態(プロトネーション状態)の標準化はこの機能の対象外。呼び出し側が
[`molstd.standardize_smiles`](../src/molstd)等で事前に行う前提。

現時点は最小限のスコープ(1配座のみ、コンフォマー探索・複数プロトネーション状態の列挙なし)で実装している。

## CLI仕様

なし。`pf`コマンド化は現在保留中([README.md](../README.md)参照)、Pythonから直接呼び出す。

```python
from ligandprep import prepare_ligand_pdbqt

output_path = prepare_ligand_pdbqt(smiles="CC(=O)Nc1ccc(O)cc1", name="acetaminophen", output_path="data/ligands/acetaminophen.pdbqt")
```

### 処理内容(`src/ligandprep/embed.py` の `prepare_ligand_pdbqt()`)

1. `Chem.MolFromSmiles(smiles)` でパースする。パース失敗は`ValueError`。
2. `Chem.AddHs(mol)` で明示的水素を付加する。
3. `AllChem.EmbedMolecule(mol, ETKDGv3(randomSeed=0xF00D))` で3D配座を1つ生成する(再現性のため固定シード)。
   生成失敗は`ValueError`。
4. `AllChem.MMFFOptimizeMolecule(mol)` でMMFF94により構造最適化する。
5. meeko `MoleculePreparation().prepare(mol)` でGasteiger電荷・原子タイプを割り当てる
   (`MoleculeSetup`のリストを返すが、単一配座入力なので`[0]`のみ使う)。
6. `PDBQTWriterLegacy.write_string(molsetup)` でPDBQT文字列に変換し`output_path`へ書き出す。変換失敗は`ValueError`。

## 実装ファイル

- `src/ligandprep/embed.py` — SMILES→3D配座→PDBQT変換(ligandprep固有のドメインロジック)

## 依存パッケージ

- `rdkit>=2026.03.5`(`pharmoforge`環境の中核ライブラリ、既存)。
- `meeko`(PDBQT変換)。conda-forgeのビルドがpython<3.14までにしか対応していないため`pharmoforge`環境に
  pipでインストールする(PyPIのwheelはpure Pythonでpython 3.14でも動作を確認済み)。依存の`gemmi`も同様の
  理由でpipインストール(PyPIにcp314向けmanylinux wheelあり)。詳細は[README.md](../README.md#依存パッケージ)参照。

## テスト

ネットワークアクセスを伴わないため、実データ(小分子SMILES)を使った実処理での検証を基本とする。

```bash
pytest tests/ligandprep
```

## 動作例

`docking`機能([cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション7)で、
ChEMBLから収集した化合物SMILESをドッキング直前にPDBQT化するために使っている。

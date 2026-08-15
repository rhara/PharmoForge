# proteinprep

PDBやAlphaFold DB由来の蛋白構造を、ドッキング・MDに使える最低限のクオリティーまで修復する機能。
現時点では最小限のスコープ(欠損原子の補完・プロトン化)のみを実装しており、
HETATM分類(リガンド抽出等)やAlphaFold構造のpLDDTに基づく末端トリムは今後の拡張予定。

## 使い方

```bash
pf prep-protein <入力構造ファイル> --output <出力PDBファイル> [--mode dock|md] [--ph <pH>]
```

```bash
pf prep-protein data/9csk.cif --output data/9csk_dock.pdb --mode dock
pf prep-protein data/P61626.cif --output data/P61626_md.pdb --mode md --ph 7.4
```

入力構造は[`pf fetch structure=...`](../fetcher/README.md)(RCSB PDB)や
[`pf fetch structure-af=...`](../fetcher/README.md)(AlphaFold DB)の出力をそのまま使える。

### 処理内容

1. PDBFixerで欠損残基・欠損原子を検出し補完する。
2. `--mode`(既定`dock`)で水素原子の扱いを切り替える。
   - `dock`: 水素原子を付加しない。
   - `md`: `--ph`(既定7.0)で指定したpHでプロトン化し、水素原子を付加する。
3. リガンドや水などのヘテロ原子は変更・除去しない(分類・整形は今後の拡張)。
4. 出力は常にPDB形式。

## 実装方針

- 構造の修復は[OpenMM](https://openmm.org/)付属の[PDBFixer](https://github.com/openmm/pdbfixer)を利用する(`src/proteinprep/repair.py`)。
- 構造ファイルの取得(RCSB PDB/AlphaFold DB)は本機能の対象外。共通パッケージ[`src/rcsb`](../src/rcsb)・[`src/afdb`](../src/afdb)を使う[`fetcher`](../fetcher/README.md)機能を利用する。

## テスト

```bash
pytest tests/proteinprep
```

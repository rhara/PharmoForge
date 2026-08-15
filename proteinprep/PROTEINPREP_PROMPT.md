# proteinprep 実装記録

このドキュメントは `proteinprep` 機能を再現するための仕様記録。

## 目的

PDB(RCSB)またはAlphaFold DB由来の蛋白構造を、ドッキング・MDに使える最低限のクオリティーまで修復する。
PharmoForgeにおける「蛋白構造の準備」機能を担う。

現時点は最小限のスコープで実装している。以下は今後の拡張候補で未実装:

- HETATM分類(water/additive/cofactor/リガンド自動抽出→SDF)
- AlphaFold DB構造のpLDDTに基づく低信頼度末端領域のトリム
- 鎖ごとの残基番号振り直し・TER整形・ジスルフィド(CYX)リネーム

## CLI仕様

```
pf prep-protein <入力構造ファイル> --output <出力PDBファイル> [--mode dock|md] [--ph <pH>]
```

- `<入力構造ファイル>` はPDBFixerが読める形式(PDB/CIF)。存在しないパスはエラー。
- `--output` / `-o` は必須。
- `--mode`(既定`dock`): `dock` は水素原子を付加しない。`md` は `--ph` で指定したpHでプロトン化する。
- `--ph`(既定`7.0`): `--mode md` の場合のみ有効。

### 処理内容(`src/proteinprep/repair.py` の `repair_structure()`)

1. `PDBFixer(filename=...)` で入力構造を読み込む。
2. `findMissingResidues()` → `findMissingAtoms()` → `addMissingAtoms()` で欠損残基・欠損原子を補完する。
   非標準残基の置換(`findNonstandardResidues`/`replaceNonstandardResidues`)やヘテロ原子の除去
   (`removeHeterogens`)は行わない(リガンド等を意図せず改変・除去しないため)。
3. `mode=md` の場合のみ `addMissingHydrogens(ph)` で指定pHのプロトン化を行う。`mode=dock` では水素原子を付加しない。
4. `openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, ..., keepIds=True)` でPDB形式に書き出す。

## 実装ファイル

- `src/proteinprep/repair.py` — PDBFixerによる欠損原子補完・プロトン化(proteinprep固有のドメインロジック)
- `src/proteinprep/cli.py` — `pf prep-protein` サブコマンド

構造データ取得(RCSB PDB/AlphaFold DB)は本機能の対象外。[`fetcher`](../fetcher/FETCHER_PROMPT.md)機能
(共通パッケージ`src/rcsb`・`src/afdb`)を利用する。

## 依存パッケージ

- `openmm` / `pdbfixer`(conda-forge)。`rdkit>=2026.03.5`との互換のため、まずは共通環境`pharmoforge`(Python 3.14)への
  インストールを試み、conda-forgeにビルドが存在することを確認済み(openmm 8.5.2 py314ビルドあり)。

## テスト

ネットワークアクセスを伴わないため、実データ(最小限の合成PDB)を使った実処理での検証を基本とする。

```bash
pytest tests/proteinprep
```

## 動作例(サンプルデータ)

ヒトリゾチームC(UniProt: P61626)のAlphaFold DB予測構造を題材にした一連の流れ:

```bash
pf fetch structure-af=P61626 --type=cif --output data/P61626.cif
pf prep-protein data/P61626.cif --output data/P61626_dock.pdb --mode dock
pf prep-protein data/P61626.cif --output data/P61626_md.pdb --mode md --ph 7.4
```

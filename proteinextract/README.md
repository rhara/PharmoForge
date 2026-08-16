# proteinextract

構造ファイル(PDB/CIF)から指定チェーンを抽出し、必要に応じて水分子を除去する機能。
ドッキングやMD等の下流処理向けに、複合体構造から必要な部分だけを切り出す用途を想定。

## 使い方

```bash
pf protein-extract <入力構造ファイル> [--chains=<チェーンID(カンマ区切り)>] [--remove-water] --output <出力構造ファイル>
```

```bash
pf protein-extract data/cdk2/2CCH.cif --chains=A --remove-water --output data/cdk2/2CCH_main.cif
pf protein-extract data/tyk2/6NZP.cif --remove-water --output data/tyk2/6NZP_nowater.pdb
```

- `--chains`(省略時は全チェーン): 抽出するチェーンIDをカンマ区切りで指定する。
- `--remove-water`(既定OFF): 水分子(`HOH`等)を除去する。
- `--output` / `-o`: 出力先。入力・出力とも拡張子(`.pdb`/`.cif`)で形式を自動判別するため、
  PDB↔CIF間の変換にも使える。
- 指定したチェーン・条件で1原子も選択されない場合は`ValueError`。

CIFファイルのチェーンIDは[`auth_asym_id`](https://mmcif.wwpdb.org/)(PyMOLやRCSBのWebサイトで
実際に表示されるチェーンID)を使う([`structio.parse_structure`](../API.md#srcstructio)が
`unite_chains=True`で読み込む)。ProDyの既定(`label_asym_id`)は同じauthチェーンに属する
水分子・リガンド等を別チェーンに細分化するため、`--chains`にPyMOL等で見えるIDをそのまま
指定できるようにするための対応。またCIF出力時、ProDyの制約によりチェーンIDが入力時と異なる
ラベルに振り直されることがある(原子の対応関係・座標自体は変化しない)。

入力構造は[`pf fetch structure=...`](../fetcher/README.md)(RCSB PDB)の出力等をそのまま使える。

## 実装方針

- `src/proteinextract/extract.py`の`extract_structure()`が本体。構造の読み書きは共通パッケージ
  [`src/structio`](../API.md#srcstructio)(ProDy)を使う。
- チェーン抽出・水除去はいずれもProDyの選択式(`chain <ID...>`、`not water`)で行う
  (`--chains`省略時は`all`)。

## テスト

```bash
pytest tests/proteinextract tests/structio
```

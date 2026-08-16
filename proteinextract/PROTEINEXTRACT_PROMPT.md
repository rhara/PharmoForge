# proteinextract 実装記録

このドキュメントは `proteinextract` 機能を再現するための仕様記録。

## 目的

構造ファイル(PDB/CIF)から指定チェーンを抽出し、必要に応じて水分子を除去して保存する。
ドッキング・MD等の下流処理向けに、複合体構造から必要な部分だけを切り出す。

## CLI仕様

```
pf protein-extract <入力構造ファイル> [--chains=<チェーンID(カンマ区切り)>] [--remove-water] --output <出力構造ファイル>
```

- `<入力構造ファイル>` はProDyが読めるPDB/CIF。存在しないパスはエラー。
- `--chains`(省略時は全チェーン): カンマ区切りのチェーンIDリスト。
- `--remove-water`(フラグ、既定OFF): 水分子を除去する。
- `--output` / `-o` は必須。拡張子(`.pdb`/`.cif`)で出力形式を判別。

### 処理内容(`src/proteinextract/extract.py`の`extract_structure()`)

1. `structio.parse_structure()`で入力構造を読み込む(拡張子で自動判別)。
2. 選択式を組み立てる: `chains`指定時は`"chain " + " ".join(chains)`、`remove_water=True`時は
   `"not water"`を追加し、`" and "`で連結する(両方省略時は`"all"`)。
3. `AtomGroup.select(selection)`で抽出する。1原子も選択されない場合は`ValueError`。
4. `structio.write_structure()`で出力する(拡張子で自動判別)。

## 実装ファイル

- `src/proteinextract/extract.py` — 抽出ロジック(proteinextract固有のドメインロジック)
- `src/proteinextract/cli.py` — `pf protein-extract`サブコマンド
- `src/structio/io.py` — 構造ファイルの読み書き(アトミックパッケージ、詳細は
  [API.md](../API.md#srcstructio))

## 動作検証で判明した知見: CIFのチェーンID(`label_asym_id` vs `auth_asym_id`)

ユーザーの依頼例(`pf protein-extract data/cdk2/2CCH.cif --chains=A,G,H,I,J ...`)で実データ
検証したところ、`2CCH.cif`は`label_asym_id`基準で17種のチェーン(`A`〜`Q`)を持つのに対し、
`unite_chains=True`で`auth_asym_id`基準にすると6種(`A`〜`F`)しかなく、依頼例の`G`〜`J`は
`auth_asym_id`には存在しないことが判明。当初はユーザーが(PyMOL等での表示と同じ)`label_asym_id`
基準のチェーンIDを想定していると解釈し、ProDyの`parseMMCIF()`の既定(`unite_chains=False`、
`label_asym_id`をチェーンIDとして扱う)をそのまま使う設計にしていた。

また、`--chains`と`--remove-water`を組み合わせて`writeMMCIF()`で書き出すと、出力ファイルを
再度読み込んだ際のチェーンIDが元のラベル(`A,G,H,I,J`)から振り直されたラベル(`A,B,C,D,E,F,G`)に
変わることを確認した(ProDyの`writeMMCIF()`はBiopythonの`MMCIFIO`経由で書き出しており、
`label_asym_id`を内部的に再割り当てするとみられる)。原子の対応関係・座標・原子数自体は正しいが、
CIF出力時にチェーンIDのラベルは保証されない。この点はユーザーに確認し、「出力の際、チェーンIDが
変わっても構わない」との回答を得たため、追加の対処はせずそのままとした。PDB出力
(`writePDB()`)では元のチェーンIDが保持されることを確認済み。

### 方針転換: `label_asym_id` → `auth_asym_id`(`unite_chains=True`)

後日、ユーザーから改めて「`pf protein-extract`のchain idの出力がおかしい、authのチェーンに
合わせてほしい」との要望を受けた。`structio.parse_structure()`は`proteinextract`だけでなく
`alignview`/`structfit`/`sequencealign`(`seqextract`経由)からも共通で使われるアトミックな
読み込み関数であるため、`proteinextract`だけでなく**この関数自体**を`unite_chains=True`に
変更する方針をユーザーに確認の上で決定(全コマンドで一貫してauthチェーンIDになる)。

これにより上記の「`--chains`にPyMOL等で見えるIDをそのまま指定できない」問題(`2CCH.cif`で
リガンド由来の`label_asym_id`(`G`〜`Q`)が本来の著者チェーン`A`〜`F`と食い違う点)を解消した。
この変更に伴い、`label_asym_id`基準の古いチェーン選択例が残っていたREADME・本ドキュメント・
`src/proteinextract/cli.py`のdocstring例(`--chains=A,G,H,I,J`)を、実データで確認した
`auth_asym_id`(`--chains=A`)に合わせて更新した(下記「動作例」参照)。

`my_examples.sh`はユーザー自身が管理する作業用スクリプトのため、Claudeからは読み込み・編集
しない方針(ユーザーの指示)。同ファイルに残る`label_asym_id`基準のチェーン選択例
(2CCHの`--chains=A,G`等)は、実行時にauth基準の`--chains=A`相当に自動的に解釈される
(存在しない`G`等は単に無視され、`A`のみが選択される。古い出力ファイル名との対応がずれる
可能性はあるが、選択結果自体は壊れない)。

## テスト

ネットワークアクセスを伴わないため、合成PDBによる実処理での検証を基本とする。

```bash
pytest tests/proteinextract tests/structio
```

## 動作例(サンプルデータ)

CDK2の結晶構造(複合体)から、蛋白(authチェーンA。旧label基準のG〜J相当のリガンドも
auth基準ではAに含まれる)のみを水を除いて抽出する:

```bash
pf protein-extract data/cdk2/2CCH.cif --chains=A --remove-water --output data/cdk2/2CCH_main.cif
```

# alignview 実装記録

このドキュメントは `alignview` 機能を再現するための仕様記録。

## 目的

複数のPDB/CIF構造をPyMOLで開き、先頭に指定した構造を基準に他の構造をアラインメント(重ね合わせ)して
表示する。AlphaFold予測構造と複数の結晶構造を素早く見比べたい、といった構造比較のためのビューア機能。

## CLI仕様

```
pf align-view <構造ファイル1> <構造ファイル2> [<構造ファイル3> ...] [--method align|super|cealign|number] [--align-margin <残基数>] [--pymol-env <env名>]
```

- 構造ファイルは1つ以上必須。存在しないパスはエラー。
- 1つめが基準(target)、2つめ以降がそれぞれ基準にアラインされる(mobile)。1つだけの場合はアラインせず開く。
- `--method`(既定`align`): `align`|`super`|`cealign`|`number`。
- `--align-margin`(既定`20`、`--method align`時のみ有効): 下記の範囲絞り込みのマージン(残基数)。
- `--pymol-env`(既定`pymol`): PyMOLがインストールされたconda/mamba環境名。

### 処理内容(`src/alignview/view.py`)

1. `_unique_object_names()`: 各ファイルのstem(拡張子なし)をPyMOLオブジェクト名にする。重複時は
   `_2`,`_3`...を付与して一意化する。
2. `build_pymol_script()`: 以下の内容の`.pml`スクリプト本文を組み立てる。
   1. 各構造を`load <絶対パス>, <オブジェクト名>`で読み込む。
   2. `show cartoon`でカートゥーン表示を明示する。
   3. オブジェクトごとに固定パレット(green/cyan/magenta/yellow/salmon/skyblue/orange/purple、
      8個超は繰り返し)で`color`する。
   4. 先頭オブジェクトを基準に、2つめ以降を`--method`に応じたコマンドでアラインする。
      - `super`: `super <mobile>, <target>`
      - `cealign`: `cealign <target>, <mobile>`(引数順が逆)
      - `align`: 単純な`align <mobile>, <target>`ではなく、`_RANGE_RESTRICTED_ALIGN_TEMPLATE`による
        `python ... python end`埋め込みブロックを生成する(詳細は下記「`align`の範囲絞り込み」)。
      - `number`: PyMOL起動前(スクリプト組み立て時)に`structfit.fit_by_residue_number()`を呼び出し、
        実際の構造ファイルから残基番号ベースの剛体変換(4x4行列)を計算し、`_NUMBER_TRANSFORM_TEMPLATE`
        による`python ... python end`ブロックで`cmd.transform_object(mobile, matrix, transpose=0)`
        として適用する(詳細は下記「`number`メソッド」)。計算に失敗した場合(共通の残基番号なし等)は
        その構造のアラインをスキップし、警告を出力する(全体は継続)。
   5. `zoom all`で全体が収まるようにビューを調整する。
3. `launch_alignment_view()`: スクリプトを`tempfile`に書き出し、
   `mamba run -n <pymol_env> pymol <script>`をサブプロセスとして起動する(GUIモード、
   ウィンドウを閉じるまでブロック)。実行後、一時ファイルを削除する。`mamba`/`conda`が
   PATH上に見つからない場合は`RuntimeError`。

## 実装ファイル

- `src/alignview/view.py` — PyMOLスクリプト組み立て・起動(alignview固有のドメインロジック)
- `src/alignview/cli.py` — `pf align-view`サブコマンド
- `src/structfit/fit.py` — 残基番号ベースの剛体重ね合わせ計算(アトミックパッケージ、詳細は
  [API.md](../API.md#srcstructfit))

## 依存パッケージ・実行環境

- `pharmoforge`環境自体は`pymol`に依存しない(外部プロセスとして`mamba run`経由で呼ぶだけ)。
- PyMOL(`pymol-open-source`、conda-forge)は`rdkit>=2026.03.5`ピンとの競合により専用環境
  (既定名`pymol`)にインストールする。
  ```bash
  mamba create -n pymol -c conda-forge pymol-open-source
  ```
- ProDy(`--method number`が使用)はconda-forgeでのビルド提供が不安定なため、`pharmoforge`環境に
  `pip`でインストールする([pyproject.toml](../pyproject.toml)の依存関係に記載)。

## 動作検証で判明した知見: `align`の範囲絞り込み

実装時、実データ(下記サンプル)で単純な`align <mobile>, <target>`(配列ベース、PyMOL標準)を
試したところ、一部構造でRMSDが13.9Å・26.4Åという明らかな失敗(対応残基の取り違え)が発生した。
TYK2はJH1(キナーゼ)/JH2(偽キナーゼ)という相同性の高い2ドメインを持ち、JH2のみをカバーする
結晶構造(6NZP・5C03)をAlphaFold全長モデルに配列アラインメントすると、配列類似度の高いJH1側に
誤対応することがあると判明。

`--method cealign`(配列非依存の構造ベース)ではこの誤対応は起きず全構造でRMSD 1.06〜1.45Åと
妥当な結果になったが、ユーザーからは「同じ蛋白なのでcealignではなくalignを使いたい」との要望が
あった。そこで、**同一蛋白の構造間ではPDBの残基番号が(UniProt基準等で)揃っている**という前提を
利用し、`align`実行前にtargetの探索対象を`resi <mobileの残基番号レンジ±margin>`に絞り込む方式を
採用した(`_RANGE_RESTRICTED_ALIGN_TEMPLATE`、PyMOLの`python`ブロックとして埋め込む)。これにより
配列アラインメントの候補が他ドメインに及ばなくなり、誤対応が解消される。同じ入力でRMSD
0.541Å・0.654Å・0.547Åと、絞り込みなしの`cealign`より高精度な結果になることを確認した
(配列対応が正確に取れる場合、構造ベースより配列ベースの方が高精度なフィットになりやすいため)。

絞り込みで対応原子が1つも得られない場合(残基番号体系が構造間で揃っていない、異なる蛋白同士の
比較等)は、`python`ブロック内で例外を捕捉し、絞り込みなしの`align`に自動フォールバックする。
残基番号が揃っていないことが最初から分かっている場合は`--method cealign`を使う。

既定値は`align`(この範囲絞り込み込み)のまま、同一蛋白の複数構造比較という主用途に最適化した。

## 動作検証で判明した知見: `number`メソッドとPyMOL自前実装の限界

ユーザーから「残基番号を明示的に使う4つめの手法として`--method number`を追加したい」との要望があった。
当初はPyMOLの`python`ブロック内で完結させる実装(`cmd.get_chains()`で鎖を選び、`cmd.iterate`で
残基番号集合を求め、`resi <番号リスト>`選択+`cmd.pair_fit()`で重ね合わせ)を試みたが、以下の問題に
遭遇し断念した。

1. `cmd.get_chains()`の返す順序が実行ごとに安定せず、`[0]`で選んだ鎖が(欠損の多い方の)NCSコピーに
   なることがあった(`sorted()`で決定的に選ぶよう修正しても解決しなかった別問題)。
2. 対応する残基数が一致していても、`resi <番号リスト>`選択で得られる原子の内部順序が両オブジェクト間で
   一致する保証がなく、`cmd.pair_fit()`が誤った対応でRMSD 20Å超という失敗をすることがあった
   (mmCIF由来のオブジェクトでは原子の格納順序が`resi`の昇順と一致しない場合があるとみられる)。
3. 誤対応を避けるため残基ごとに1原子ずつ明示ペアを`cmd.pair_fit()`に渡す(数百残基分の引数)方式を
   試したところ、PyMOL(pymol-open-source)がセグメンテーション違反で異常終了した。

ユーザーから「必要なときはProDyまたはBioPythonを使ってください」との指示を受け、方針を転換。
構造のパース・残基対応付け・重ね合わせ計算そのものは、PyMOLではなく確立された構造生物学ライブラリで
`pharmoforge`環境側(PyMOL起動前)に行い、その結果(4x4剛体変換行列)だけをPyMOLに渡して
`cmd.transform_object()`で適用する設計に変更した。ライブラリはBioPython(`Bio.PDB`)とProDyの両方で
試作し、結果は同等だったが、ユーザーの指定により最終的にProDy(`pip install prody`、
[pyproject.toml](../pyproject.toml)に追記)を採用した。この経緯から、新規パッケージ
[`src/structfit`](../API.md#srcstructfit)を切り出した(構造間の残基番号ベース重ね合わせという、
alignview以外でも再利用しうるアトミックな処理のため)。

`cmd.transform_object(name, matrix, transpose=0)`が期待する`matrix`の規約(4x4、行優先、
`v' = matrix @ [x, y, z, 1]`)は、実際に変換を適用して特定残基のCA座標を検証し確認した
(構造上の1残基で変換後座標とtarget側座標の距離が0.43Åとなり、全体RMSD(0.79Å)と整合することを確認)。

同じ入力データ(TYK2)で`--method number`を実行するとRMSD 0.79〜0.93Åとなり、`--method align`
(0.54〜0.65Å、外れ値除去による精密化あり)よりはやや大きいが、`cealign`(1.06〜1.45Å)より高精度で、
配列アラインメントを一切経由しない最も直接的な対応付けとして機能することを確認した。

## テスト

ネットワークアクセスを伴わないため、`subprocess.run`(および`--method number`では
`structfit.fit_by_residue_number`)をモックした単体テストを基本とする。`structfit`自体は
合成PDBによる実処理での検証を行う([tests/structfit](../tests/structfit))。

```bash
pytest tests/alignview tests/structfit
```

## 動作例(サンプルデータ)

TYK2のAlphaFold DB予測構造と複数のPDB結晶構造を重ね合わせて表示する:

```bash
pf fetch structure=TYK2_HUMAN --af --type=cif --output data/tyk2/TYK2_HUMAN_af.cif
pf fetch structure=6NZP --output data/tyk2/6NZP.cif
pf fetch structure=4OLI --output data/tyk2/4OLI.cif
pf fetch structure=5C03 --output data/tyk2/5C03.cif
pf align-view data/tyk2/TYK2_HUMAN_af.cif data/tyk2/6NZP.cif data/tyk2/4OLI.cif data/tyk2/5C03.cif
pf align-view data/tyk2/TYK2_HUMAN_af.cif data/tyk2/6NZP.cif --method number
```

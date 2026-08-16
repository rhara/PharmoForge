# alignview

複数のPDB/CIF構造をPyMOLで開き、先頭に指定した構造を基準に他の構造をアラインメント(重ね合わせ)して
表示する機能。構造比較(AlphaFold予測構造と複数の結晶構造の見比べ等)を素早く行うためのビューア。

## 使い方

```bash
pf align-view <構造ファイル1> <構造ファイル2> [<構造ファイル3> ...] [--method align|super|cealign|number] [--align-margin <残基数>] [--pymol-env <env名>]
```

```bash
pf align-view data/tyk2/TYK2_HUMAN_af.cif data/tyk2/6NZP.cif data/tyk2/4OLI.cif data/tyk2/5C03.cif
pf align-view data/tyk2/TYK2_HUMAN_af.cif data/tyk2/6NZP.cif --method number
pf align-view data/9csk.cif data/1abc.cif --method cealign
pf align-view --indir data/cyp P08604_AF 1PQ2_ad 3IBD_abcde 3NXU_abh
pf align-view --indir data/cyp P08604_AF 1PQ2_ad --indir data/other 9XYZ
```

- 1つめの構造ファイルが基準(target)、2つめ以降がそれぞれ基準にアラインされる(mobile)。
- 構造ファイルが1つだけの場合はアラインは行わず、そのままPyMOLで開く。
- 各構造はファイル名(拡張子なし)をオブジェクト名としてPyMOLに読み込まれ、見分けやすいよう自動で配色される。
- `--indir DIR`: 繰り返し指定可能。以降のファイル名(拡張子省略可、`.cif`優先、次に`.pdb`)を
  `DIR`配下から解決する。`/`を含む指定(または絶対パス)は`--indir`によらずカレントディレクトリ
  相対 or 絶対パスとして扱う([`structio.resolve`](../API.md#srcstructio)、`sequence-align`と共用)。
- `--method`(既定`align`): アラインメント手法。同一蛋白の構造間ではPDBの残基番号が(UniProt基準等で)
  揃っている前提のもと、`align`/`number`はいずれも配列アラインメントではなく残基番号の対応付けを利用する。
  - `align`: targetの探索範囲をmobileの残基番号レンジ(`--align-margin`分の余裕を加えた範囲)に
    絞り込んだ上でPyMOL標準のalign(配列ベース)を実行する。これにより、相同性の高い別ドメイン
    (例: キナーゼ/偽キナーゼの対)に誤対応することを防ぐ。絞り込みで対応が取れない場合は自動的に
    絞り込みなしのalignにフォールバックする。外れ値を反復的に除外して精密化するため、対応する
    範囲では最も高精度になりやすい。
  - `number`: 配列アラインメントを一切行わず、残基番号が一致するCA原子同士を直接対応付けて
    剛体重ね合わせを行う([`structfit`](../API.md#srcstructfit)/ProDyによる計算)。最も直接的・
    決定的な対応付け。複数鎖を含む構造では、共通残基番号数が最大になる鎖の組を自動選択する。
  - `super`/`cealign`: 配列非依存の構造ベースアラインメント。異なる蛋白同士の比較や、構造間で
    残基番号が揃っていない場合に使う。
- `--align-margin`(既定`20`、`--method align`時のみ有効): targetの探索範囲をmobileの残基番号
  レンジ±この値に絞り込む際のマージン。
- `--pymol-env`(既定`pymol`): PyMOLをインストールした専用conda/mamba環境名(下記参照)。

入力構造は[`pf fetch structure=...`](../fetcher/README.md)(RCSB PDB)や
[`pf fetch structure=... --af`](../fetcher/README.md)(AlphaFold DB)の出力、または
[`pf prep-protein`](../proteinprep/README.md)の出力をそのまま使える。

## PyMOLの実行環境

PyMOL(`pymol-open-source`)はrdkitのバージョン要件(`rdkit>=2026.03.5`)が共通環境`pharmoforge`と
競合するため、専用のconda/mamba環境(既定名`pymol`)にインストールする前提。`pf align-view`は
`mamba run -n <env名> pymol <生成したスクリプト>`という形で外部プロセスとして起動する
(`pharmoforge`環境自体は`pymol`パッケージに依存しない)。

```bash
mamba create -n pymol -c conda-forge pymol-open-source
```

`--method number`は[`structfit`](../API.md#srcstructfit)(ProDy)を使い、PyMOL起動前に
`pharmoforge`環境側で重ね合わせを計算する。ProDyは`pharmoforge`環境に`pip`でインストール済み
([pyproject.toml](../pyproject.toml)参照)。

## 実装方針

- `src/alignview/view.py`の`build_pymol_script()`で、構造の読み込み・配色・アラインメントを
  行うPyMOLスクリプト(`.pml`)を組み立てる。`--method number`の場合はここで(PyMOL起動前に)
  `structfit.fit_by_residue_number()`を呼び出し、実際に構造ファイルを読んで剛体変換を計算する。
- `launch_alignment_view()`がスクリプトを一時ファイルに書き出し、`mamba run -n <env> pymol <script>`で
  PyMOLをGUIモードで起動する(処理はPyMOLウィンドウを閉じるまでブロックする)。実行後、一時ファイルは削除する。
- `--indir`解決ロジック(拡張子省略時の自動補完等)は[`structio.resolve`](../API.md#srcstructio)を
  `sequence-align`と共用する。

## テスト

```bash
pytest tests/alignview tests/structfit
```

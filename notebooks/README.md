# notebooks

`src`以下の関数(`API.md`参照)を実際に使い、テーマを持った創薬支援タスクをJupyter notebookとして
まとめる場所。`pf`コマンド(CLI)化は前提とせず、notebook自体が活用例・成果物になる。

## 環境

`pharmoforge`環境に`jupyterlab`/`ipykernel`(conda-forge)をインストール済み。

```bash
mamba install -n pharmoforge -c conda-forge jupyterlab ipykernel
python -m ipykernel install --user --name pharmoforge --display-name "Python (pharmoforge)"
```

`pharmoforge`は`pip install -e .`でエディタブルインストールされているため、notebook内でも
`from pocket import run_fpocket`のように`src`以下のパッケージをそのままimportできる
(パス操作は不要)。Notebookのカーネルには"Python (pharmoforge)"を選択する。

```bash
jupyter lab
```

## 命名・構成の方針

- 1テーマ1notebook。ファイル名は`<連番>_<テーマ>.ipynb`のように内容がわかる名前にする。
- 出力先ディレクトリ(構造ファイル・画像等の生成物)は、notebookごとに専用のサブディレクトリ
  (例: `notebooks/data/<テーマ名>/`)を使う。複数notebookで`outdir`等をベタ書きの共通名にすると
  cwdを共有した際に衝突するため、必ずnotebook固有の名前にする。

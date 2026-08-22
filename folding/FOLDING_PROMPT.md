# folding 実装記録

このドキュメントは `folding` 機能を再現するための仕様記録。関数シグネチャ・挙動の正は
[API.md](../API.md#srcfolding)を参照。

## 目的

[Boltz-2](https://github.com/jwohlwend/boltz)(MIT/NVIDIA、AlphaFold3系のオープンソース生体分子構造
予測モデル)による蛋白構造予測。PharmoForgeにおける「フォルディング予測」機能を担う。

現時点は最小限のスコープで実装している。以下は今後の拡張候補で未実装:

- `pf`コマンド化(現在は`pharmoforge`ライブラリ関数として直接呼び出す運用、[README.md](../README.md)参照)
- 蛋白-リガンド共フォールディング(`Ligand`/`predict_affinity`)、共有結合複合体、DNA-蛋白複合体
- `pocket`/`contact`/`bond`拘束(NVIDIAホスト型APIは対応済みだが未実装)
- ローカルGPU実行(下記「経緯」参照、現状ホスト型APIのみ)

## 経緯: ローカルGPUではなくNVIDIAホスト型APIを使う理由

Boltz-2の実測メモリ使用量は単純なモノマー予測でも数GB〜(文献ベンチマークではL40S 48GB環境で
構造予測に約11GB使用の報告あり)、このプロジェクトの検証機(NVIDIA GeForce GTX 1660 Ti、VRAM 6GB)
では不足すると判断し、ローカル実行は行わない方針とした(ユーザー指示)。代わりにNVIDIAのホスト型API
(`https://health.api.nvidia.com`)を使う。認証は`.env`の`NVIDIA_API_KEY`(`nvapi-...`形式)。

`.env`はリポジトリに存在するが誤ってgitignoreされていなかったため、この対応の過程で
`.gitignore`に追加した(秘密情報の誤commit防止)。

## 経緯: Pythonクライアントの選定・既知の不具合

NVIDIA公式の[`boltz2-python-client`](https://github.com/NVIDIA/digital-biology-examples/tree/main/examples/nims/boltz-2)
(PyPI、conda-forgeには無いためpipインストール)を使う。純粋なHTTPクライアント(httpx/aiohttp/pydantic等)
でGPU/CUDAへの依存が無く、`pharmoforge`環境に競合なくインストールできることを確認済み。

**既知の不具合**: `MSASearchClient`(MSA検索NIM呼び出し)の内部パス定数
(`_MONOMER_PATH = "/biology/colabfold/msa-search/predict"`等)に`/v1`プレフィックスが欠落しており、
`configure_msa_search()`のドキュメント例通り`endpoint_url="https://health.api.nvidia.com"`を渡すと
404になる(実測確認済み)。`endpoint_url="https://health.api.nvidia.com/v1"`(`/v1`を含める)ことで
回避できる。`src/folding/boltz.py`の`_MSA_BASE_URL`定数にこの回避策を実装済み。

構造予測本体(`Boltz2Client.predict`)のエンドポイントは`base_url`に`/v1`を含めない
(`https://health.api.nvidia.com`)のが正しく、こちらは不具合なし(内部で`predict_url`が
適切に`/v1/biology/mit/boltz2/predict`を組み立てる)。

## CLI仕様

なし。Pythonから直接呼び出す(使い方は[README.md](../README.md#使い方)参照)。

### 処理内容(`src/folding/boltz.py`)

1. `search_msa()`: `MSASearchClient` + `MSASearchIntegration.search_and_save()`でNVIDIA MSA Search
   NIM(colabfold)にリクエストし、a3m形式のアラインメントをファイルに保存する。
2. `predict_structure()`:
   - `msa_path`が指定されていれば、その内容を`AlignmentFileRecord`にラップし
     `{"msa_search": {"a3m": AlignmentFileRecord(...)}}`の形(Boltz-2 APIが要求するネスト構造、
     `boltz2_client`内部の`search_and_prepare_for_boltz()`と同じ形)で`Polymer.msa`に渡す。
   - `templates`が指定されていれば、各テンプレートファイルの内容を`StructuralTemplate`
     (`structure`/`format`/`chain_id`/`name`)に変換し`Polymer.structural_templates`に渡す
     (最大4件、蛋白のみ)。
   - `Boltz2Client(endpoint_type=EndpointType.NVIDIA_HOSTED).predict(request)`を呼び、
     `PredictionResponse.structures`(モデルごとのmmCIF文字列)を`<output_dir>/<name>_model<N>.cif`
     として書き出す。

## テンプレート誘導の検証

CDK20(UniProt Q8IZL9)についてDFG-in(テンプレート無し)・DFG-out(テンプレートで誘導)の複数条件で
予測し、DFGモチーフ(CDK20の`Phe146`、`kinasemotifs.find_kinase_motifs`で検出)の位置を比較した。

測定した指標(いずれも構造内の原子間距離、フレーム非依存):
- `Phe146(CZ)` <-> 触媒`Lys33(NZ)`
- `Phe146(CZ)` <-> hinge領域`resnum 65(CA)`

### 試行1: p38 MAPK(PDB 1KV1、Type II阻害剤BIRB-796結合、CMGCグループだが別サブファミリー)

| 条件 | Phe-Lys距離(5サンプル) | Phe-hinge距離(5サンプル) |
| --- | --- | --- |
| DFG-in(テンプレート無し、MSAあり) | 10.33-11.75 Å(平均11.22) | (未測定) |
| DFG-out(1KV1テンプレート、MSAあり) | 10.73-12.19 Å(平均11.32) | 6.95-7.20 Å |
| DFG-out(1KV1テンプレート、MSA無し=単一配列モード) | 10.98-11.57 Å | 7.35-7.52 Å |

### 試行2: CDK2(PDB 5A14、Type II阻害剤LQ5結合、CDK20と同じCDKサブファミリー、配列同一性~44-45%)

より近縁な(同じCDKサブファミリーの)テンプレートに差し替えて再検証した。

| 条件 | Phe-Lys距離(5サンプル) |
| --- | --- |
| DFG-in(テンプレート無し、MSAあり) | 10.33-11.75 Å(平均11.22) |
| DFG-out(5A14テンプレート、MSAあり) | 10.80-12.41 Å(平均11.49) |
| **参照: CDK2(5A14)自身の実際のDFG-out構造** | **8.91 Å** |

CDK2(5A14)自身の実測値(8.91 Å)は、CDK20のDFG-in・DFG-out(いずれのテンプレート試行も)の
予測値群(10.3〜12.4 Å)とは明確に異なる。つまりCDK20の予測群はテンプレートの有無によらず
「実際のDFG-out状態」から遠く、**テンプレートによる誘導は確認できなかった**。

(5A14はN末端側に構造未決定領域があり、chainのauth番号にギャップがある。DFGのPhe・触媒Lysの
auth番号はたまたまCDK20と同じ(146・33)だが、これは`get_chain_sequences`の`resnums`で実際の
auth番号を確認した上での値であり、`find_kinase_motifs`が返す配列内1始まり位置をそのままauth番号
とみなすと(観測残基のみの連結配列に対する位置のため)ズレることに注意。)

### 結論

2種類のテンプレート(系統的に離れたp38、近縁なCDK2)のいずれでも、MSAの有無によらず結果は
変わらなかった。**テンプレート選択の問題ではなく、NVIDIAホスト型APIのテンプレート機構自体が
原因である可能性が高い**と考えられる:

- NVIDIAホスト型APIのテンプレート機構(`structural_templates`)にはOSS版`boltz`の
  `force`/`threshold`(RMSD拘束ポテンシャル)が無く、モデルのテンプレート注意機構による
  ソフトな誘導のみである(上記「経緯」参照)。この誘導力が、大きな構造変化(DFGモチーフの
  フリップ)を起こすには不十分だった可能性が高い。

今後の選択肢(未実施、ユーザーとの相談が必要):
- OSS版`boltz`をローカルGPUで実行し`force`/`threshold`拘束を使う(VRAM不足のため要検討)
- `recycling_steps`/`sampling_steps`を増やす、複数テンプレートを組み合わせる等のパラメータ探索

## 依存パッケージ

- `boltz2-python-client`(pip、conda-forgeに存在しない)。
- `python-dotenv`(conda-forge)。`.env`からの認証情報読み込みに使用。

## テスト

ネットワークアクセス(有償API呼び出し)を伴うため、`boltz2_client.Boltz2Client`/
`MSASearchIntegration`をモックして検証する(`docking.vina`のテストパターンを踏襲)。

```bash
pytest tests/folding
```

## 動作例(サンプルデータ)

[cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション8: CDK20(Q8IZL9)の
MSA検索、DFG-in予測、DFG-out誘導予測(p38 MAPK 1KV1テンプレート)、DFGモチーフ位置の定量比較。

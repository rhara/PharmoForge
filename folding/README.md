# folding

[Boltz-2](https://github.com/jwohlwend/boltz)による蛋白構造予測を行う機能。ローカルGPU(このプロジェクトの
検証機はVRAM不足)ではなく、NVIDIAのホスト型API(`https://health.api.nvidia.com`)を使う。

`pf`コマンド化は現在保留中(トップ[README.md](../README.md)参照)。現時点ではPythonから直接呼び出す。
関数シグネチャの正は[API.md](../API.md#srcfolding)を参照(ここでは典型的な使い方の流れのみ示す)。
実データでの使用例は[cdk20_investigation.ipynb](../notebooks/cdk20_investigation.ipynb)セクション8。

## 認証

環境変数`NVIDIA_API_KEY`(`nvapi-...`形式、[build.nvidia.com](https://build.nvidia.com)発行のAPIキー)。
PharmoForgeでは`.env`(gitignore済み)にキーを置き、`python-dotenv`で読み込む運用:

```python
from dotenv import load_dotenv
load_dotenv()  # .envのNVIDIA_API_KEYを環境変数に読み込む
```

## 使い方

```python
from folding import search_msa, predict_structure, StructureTemplate

# 1. MSA検索(1回実行し、複数の予測(DFG-in/DFG-out等)で使い回す)
msa_path = search_msa("MDQYCILGRIG...", "data/cdk20/boltz/cdk20.a3m")

# 2. テンプレートなしで予測(通常の予測、例: DFG-in)
result = predict_structure(
    "MDQYCILGRIG...", "data/cdk20/boltz/dfg_in", "cdk20_dfg_in",
    msa_path=msa_path,
)

# 3. テンプレート構造で誘導した予測(例: 別キナーゼのDFG-out構造をテンプレートに)
result_out = predict_structure(
    "MDQYCILGRIG...", "data/cdk20/boltz/dfg_out", "cdk20_dfg_out",
    msa_path=msa_path,
    templates=[StructureTemplate(structure_path="data/1KV1.cif", chain_id="A")],
)

print(result.structure_paths, result.confidence_scores)
```

### 処理内容

1. **MSA検索**(`search_msa`): NVIDIAのMSA Search NIM(colabfold)で配列に対する多重配列アラインメントを
   検索し、a3m形式でファイルに保存する。
2. **構造予測**(`predict_structure`): Boltz-2 NIMで構造を予測する。`msa_path`省略時は単一配列モード
   (精度低下、非推奨)。`templates`(最大4件、蛋白のみ)を指定すると、構造予測をそのテンプレートに
   近づけるよう誘導する。生成された各モデル(`diffusion_samples`件)をCIFとして保存する。

## テンプレートによる構造誘導の限界(実測)

CDK20のDFG-out構造を、系統的に離れたp38 MAPK(PDB 1KV1)・近縁なCDK2(PDB 5A14、同じCDK
サブファミリー、Type II阻害剤LQ5結合のDFG-out構造)の2種類のテンプレートで誘導しようと試みたが、
**いずれも実測ではDFGモチーフの構造がテンプレート無しの場合とほぼ変化しなかった**(MSA有無を
問わず)。CDK2(5A14)自身の実際のDFG-out構造を基準値として測ると、CDK20の予測群(DFG-in・
DFG-out試行とも)とは明確に異なる値であり、CDK20の予測は一貫してDFG-in寄りのままだった。
テンプレート選択の問題ではなく、NVIDIAホスト型APIのテンプレート機構自体(OSS版`boltz`CLIにある
`force`/`threshold`(RMSDを一定範囲内に強制する専用ポテンシャル)を持たない、モデルのテンプレート
注意機構によるソフトな誘導のみ)が、DFGモチーフのフリップのような大きな構造変化を起こすには
不十分である可能性が高い。詳細・実測データは[FOLDING_PROMPT.md](FOLDING_PROMPT.md#テンプレート誘導の検証)参照。

## テスト

```bash
pytest tests/folding
```

"""NVIDIA Boltz-2 NIM(ホスト型API)による蛋白構造予測。

ローカルGPU(VRAM不足)ではなくNVIDIAのホスト型API(`https://health.api.nvidia.com`)を使う。
認証は環境変数`NVIDIA_API_KEY`(`nvapi-...`形式、NVIDIA Build/NIMのAPIキー)。
"""

import asyncio
import os
from dataclasses import dataclass
from pathlib import Path

from boltz2_client import Boltz2Client, EndpointType
from boltz2_client.models import AlignmentFileRecord, Polymer, PredictionRequest, StructuralTemplate
from boltz2_client.msa_search import MSASearchClient, MSASearchIntegration

from core.logging_utils import get_logger

logger = get_logger(__name__)

_PREDICT_BASE_URL = "https://health.api.nvidia.com"
# MSA検索NIM(colabfold)の実エンドポイントは `/v1/biology/colabfold/msa-search/...` だが、
# boltz2_client.msa_search.MSASearchClientの内部パス定数(_MONOMER_PATH等)には`/v1`が
# 含まれておらず、base_urlをそのまま渡すと404になる(2025年時点のクライアントのバグと思われる)。
# base_url側に`/v1`を含めることで辻褄を合わせる。
_MSA_BASE_URL = "https://health.api.nvidia.com/v1"


def _resolve_api_key(api_key: str | None) -> str:
    key = api_key or os.environ.get("NVIDIA_API_KEY")
    if not key:
        raise ValueError(
            "NVIDIA_API_KEYが設定されていない(引数api_key、または環境変数NVIDIA_API_KEYで指定。"
            "dotenv.load_dotenv()で.envから読み込む運用を想定)"
        )
    return key


@dataclass
class StructureTemplate:
    """構造予測を誘導するテンプレート構造(1件、Boltz-2 NIMは蛋白ポリマーごとに最大4件まで)。"""

    structure_path: Path
    chain_id: str | None = None  # テンプレート構造中のどのチェーンを使うか
    name: str | None = None


@dataclass
class BoltzPrediction:
    """Boltz-2予測結果。"""

    structure_paths: list[Path]  # モデルごとのCIFファイル(信頼度降順)
    confidence_scores: list[float]
    ptm_scores: list[float]


def search_msa(
    sequence: str,
    output_path: Path,
    api_key: str | None = None,
    max_msa_sequences: int = 500,
    e_value: float = 0.0001,
) -> Path:
    """配列に対しNVIDIA MSA Search NIM(colabfold)でMSAを検索し、a3m形式で`output_path`に保存する。

    Boltz-2予測(`predict_structure`)にそのまま渡せる(`msa_path`引数)。1回の検索結果を
    DFG-in/DFG-out等の複数予測で使い回すことを想定し、検索とファイル保存を分離している。
    """
    api_key = _resolve_api_key(api_key)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    async def _run() -> Path:
        client = MSASearchClient(endpoint_url=_MSA_BASE_URL, api_key=api_key)
        integration = MSASearchIntegration(client)
        return await integration.search_and_save(
            sequence=sequence,
            output_path=output_path,
            output_format="a3m",
            max_msa_sequences=max_msa_sequences,
            e_value=e_value,
        )

    logger.info("Searching MSA for sequence (%d residues) via NVIDIA MSA Search NIM ...", len(sequence))
    result_path = asyncio.run(_run())
    logger.info("Saved MSA -> %s", result_path)
    return Path(result_path)


def predict_structure(
    sequence: str,
    output_dir: Path,
    name: str,
    msa_path: Path | None = None,
    templates: list[StructureTemplate] | None = None,
    recycling_steps: int = 3,
    sampling_steps: int = 50,
    diffusion_samples: int = 1,
    api_key: str | None = None,
    timeout: float = 1800.0,
) -> BoltzPrediction:
    """Boltz-2 NIMで蛋白構造を予測し、生成された各モデルをCIFとして`output_dir`に保存する。

    `msa_path`省略時は単一配列モード(精度が下がるため`search_msa`で事前に用意することを推奨)。
    `templates`を指定すると構造予測をそのテンプレートに近づけるようソフトに誘導する
    (OSS版`boltz`の`force`/`threshold`のような厳密な拘束ポテンシャルはNVIDIAホスト型APIの
    リクエストスキーマには存在せず、モデルのテンプレート注意機構によるソフトな誘導のみ)。
    最大4テンプレートまで。
    """
    api_key = _resolve_api_key(api_key)

    msa = None
    if msa_path is not None:
        a3m_text = Path(msa_path).read_text()
        msa = {"msa_search": {"a3m": AlignmentFileRecord(alignment=a3m_text, format="a3m", rank=0)}}

    structural_templates = None
    if templates:
        if len(templates) > 4:
            raise ValueError(f"テンプレートは最大4件まで(指定: {len(templates)}件)")
        structural_templates = [
            StructuralTemplate(
                structure=Path(t.structure_path).read_text(),
                format=(Path(t.structure_path).suffix.lstrip(".").lower() or "cif"),
                chain_id=t.chain_id,
                name=t.name,
            )
            for t in templates
        ]

    polymer = Polymer(
        id="A",
        molecule_type="protein",
        sequence=sequence,
        msa=msa,
        structural_templates=structural_templates,
    )
    request = PredictionRequest(
        polymers=[polymer],
        recycling_steps=recycling_steps,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
    )

    async def _run():
        client = Boltz2Client(
            base_url=_PREDICT_BASE_URL, api_key=api_key, endpoint_type=EndpointType.NVIDIA_HOSTED, timeout=timeout,
        )
        return await client.predict(request, save_structures=False)

    logger.info(
        "Predicting structure for %s (%d residues, %d template(s), msa=%s) via Boltz-2 NIM ...",
        name, len(sequence), len(structural_templates or []), msa_path is not None,
    )
    response = asyncio.run(_run())

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    structure_paths = []
    for i, structure in enumerate(response.structures, start=1):
        path = output_dir / f"{name}_model{i}.cif"
        path.write_text(structure.structure)
        structure_paths.append(path)

    logger.info(
        "Done: %s -> %d structure(s), confidence=%s",
        name, len(structure_paths), [round(c, 3) for c in response.confidence_scores],
    )
    return BoltzPrediction(
        structure_paths=structure_paths,
        confidence_scores=response.confidence_scores,
        ptm_scores=response.ptm_scores,
    )

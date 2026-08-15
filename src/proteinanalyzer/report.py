"""蛋白情報dictのファイル書き出し。"""

import json
from pathlib import Path

from core.logging_utils import get_logger

logger = get_logger(__name__)


def write_protein_info_json(info: dict, output: Path) -> None:
    """蛋白情報dictをJSON(インデント付き、日本語等はエスケープしない)としてoutputに書き出す。"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n")
    logger.info("Done: saved protein info (%s) to %s", info.get("accession"), output)

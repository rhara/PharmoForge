"""RCSB PDBから蛋白構造ファイルを取得する。"""

from pathlib import Path

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.{fmt}"


def fetch_structure(pdb_id: str, output: Path) -> None:
    """PDB ID(例: 9CSK)の構造ファイルをダウンロードしoutputへ保存する。

    フォーマットはoutputの拡張子(.cif / .pdb)から決める。省略時はcif。
    """
    pdb_id = pdb_id.strip().upper()
    fmt = output.suffix.lstrip(".").lower() or "cif"
    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id, fmt=fmt)

    logger.info("Downloading structure %s (%s) from RCSB PDB ...", pdb_id, fmt)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(resp.content)
    logger.info("Done: saved %s (%d bytes) to %s", pdb_id, len(resp.content), output)

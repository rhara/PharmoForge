"""RCSB PDBから蛋白構造ファイルを取得する。"""

from pathlib import Path

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.{fmt}"


def fetch_structure(pdb_id: str, output: Path, fmt: str | None = None) -> None:
    """PDB ID(例: 9CSK)の構造ファイルをダウンロードしoutputへ保存する。

    fmt(cif/pdb)を省略した場合は、outputの拡張子(.cif / .pdb)から決める。
    どちらもなければcif。
    """
    pdb_id = pdb_id.strip().upper()
    fmt = (fmt or output.suffix.lstrip(".") or "cif").lower()
    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id, fmt=fmt)

    logger.info("Downloading structure %s (%s) from RCSB PDB ...", pdb_id, fmt)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(resp.content)
    logger.info("Done: saved %s (%d bytes) to %s", pdb_id, len(resp.content), output)


def fetch_structures(pdb_ids: list[str], output_dir: Path, fmt: str) -> None:
    """複数のPDB IDをまとめてダウンロードし、`output_dir/<PDB_ID>.<fmt>` に保存する。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    total = len(pdb_ids)
    logger.info("Downloading %d structures (%s) into %s ...", total, fmt, output_dir)
    for i, pdb_id in enumerate(pdb_ids, start=1):
        pdb_id = pdb_id.strip().upper()
        logger.info("[%d/%d] %s", i, total, pdb_id)
        fetch_structure(pdb_id, output_dir / f"{pdb_id}.{fmt}", fmt=fmt)
    logger.info("Done: %d structures saved to %s", total, output_dir)

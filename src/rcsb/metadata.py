"""RCSB PDBエントリのメタデータ(実験手法・解像度等)を取得する。"""

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


def fetch_entry_info(pdb_id: str) -> dict:
    """PDB IDの実験手法・解像度等をRCSB Data APIから取得する。

    返り値: {"pdb_id", "method", "resolution"}。resolutionはX線構造以外ではNoneになりうる。
    """
    pdb_id = pdb_id.strip().upper()
    resp = requests.get(RCSB_ENTRY_URL.format(pdb_id=pdb_id), timeout=30)
    resp.raise_for_status()
    entry_info = resp.json().get("rcsb_entry_info", {})
    resolution = entry_info.get("resolution_combined")
    return {
        "pdb_id": pdb_id,
        "method": entry_info.get("experimental_method"),
        "resolution": resolution[0] if resolution else None,
    }


def fetch_entries_info(pdb_ids: list[str]) -> list[dict]:
    """複数PDB IDのメタデータをまとめて取得する。"""
    total = len(pdb_ids)
    logger.info("Fetching entry info for %d PDB entries ...", total)
    results = []
    for i, pdb_id in enumerate(pdb_ids, start=1):
        info = fetch_entry_info(pdb_id)
        logger.info(
            "  [%d/%d] %s: method=%s resolution=%s",
            i,
            total,
            info["pdb_id"],
            info["method"],
            info["resolution"],
        )
        results.append(info)
    logger.info("Done: fetched info for %d entries", total)
    return results

"""ChEMBL REST APIから活性データ(requestsベース、chembl_webresource_client不使用)を取得する。"""

import csv
import statistics
from collections import defaultdict
from pathlib import Path

import requests

from core.logging_utils import get_logger
from molstd import standardize_smiles

logger = get_logger(__name__)

CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CHEMBL_ORIGIN = "https://www.ebi.ac.uk"

AGGREGATED_FIELDS = [
    "smiles",
    "_median",
    "_mean",
    "_sd",
    "_n",
]


def fetch_activities(target_chembl_id: str, page_size: int = 1000) -> list[dict]:
    """指定したChEMBL target idについて、pChEMBL値を持つ活性データを全件取得する。"""
    url = CHEMBL_ACTIVITY_URL
    params = {
        "target_chembl_id": target_chembl_id,
        "pchembl_value__isnull": "false",
        "limit": page_size,
        "offset": 0,
    }

    records: list[dict] = []
    page = 1
    logger.info("Fetching activities for target %s (pChEMBL value required) ...", target_chembl_id)
    while url:
        resp = requests.get(url, params=params if page == 1 else None, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        activities = data.get("activities", [])
        records.extend(activities)
        logger.info("  page %d: +%d records (total %d)", page, len(activities), len(records))

        next_path = data.get("page_meta", {}).get("next")
        url = f"{CHEMBL_ORIGIN}{next_path}" if next_path else None
        page += 1

    logger.info("Done: %d activities fetched for %s", len(records), target_chembl_id)
    return records


def standardize_and_aggregate(records: list[dict]) -> list[dict]:
    """化合物構造をChEMBL Structure Pipelineに倣って標準化し、
    標準化後の構造が同じ化合物のpChEMBL値をmean/median/sdに集約する。
    """
    logger.info("Standardizing structures and aggregating pChEMBL values ...")
    groups: dict[str, list[float]] = defaultdict(list)
    skipped = 0

    for record in records:
        smiles = record.get("canonical_smiles")
        pchembl_raw = record.get("pchembl_value")
        if not smiles or pchembl_raw is None:
            skipped += 1
            continue
        try:
            pchembl_value = float(pchembl_raw)
        except (TypeError, ValueError):
            skipped += 1
            continue

        std_smiles = standardize_smiles(smiles)
        if std_smiles is None:
            skipped += 1
            continue

        groups[std_smiles].append(pchembl_value)

    if skipped:
        logger.info("  skipped %d records (missing/invalid SMILES or pChEMBL value)", skipped)
    logger.info("  %d unique standardized compounds", len(groups))

    aggregated = []
    for std_smiles, values in groups.items():
        aggregated.append(
            {
                "smiles": std_smiles,
                "_median": round(statistics.median(values), 3),
                "_mean": round(statistics.mean(values), 3),
                "_sd": round(statistics.stdev(values), 3) if len(values) >= 2 else "",
                "_n": len(values),
            }
        )
    aggregated.sort(key=lambda row: row["_median"], reverse=True)
    return aggregated


def write_activities_tsv(records: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing %d records to %s ...", len(records), output)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=AGGREGATED_FIELDS, delimiter="\t", extrasaction="ignore", restval=""
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record)
    logger.info("Done: wrote %s", output)

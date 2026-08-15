"""ChEMBL REST APIから活性データ(requestsベース、chembl_webresource_client不使用)を取得する。"""

import csv
from pathlib import Path

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CHEMBL_ORIGIN = "https://www.ebi.ac.uk"

ACTIVITY_FIELDS = [
    "molecule_chembl_id",
    "canonical_smiles",
    "target_chembl_id",
    "target_pref_name",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "assay_chembl_id",
    "assay_description",
    "document_chembl_id",
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


def write_activities_tsv(records: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing %d records to %s ...", len(records), output)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=ACTIVITY_FIELDS, delimiter="\t", extrasaction="ignore", restval=""
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record)
    logger.info("Done: wrote %s", output)

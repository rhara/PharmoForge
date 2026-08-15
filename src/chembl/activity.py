"""ChEMBL REST APIから活性データを取得する(requestsベース、chembl_webresource_client不使用)。"""

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CHEMBL_ORIGIN = "https://www.ebi.ac.uk"


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

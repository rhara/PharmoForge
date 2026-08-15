"""UniProt entry name / accession と ChEMBL target id の相互マッピング。

現時点では「識別子 -> ChEMBL target id」「識別子 -> UniProt accession」の方向のみ実装している。
逆方向の解決が必要になった時点で、この下に追加していく。
"""

import re

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
CHEMBL_TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target.json"

_CHEMBL_ID_RE = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)
_UNIPROT_ACCESSION_RE = re.compile(
    r"^([A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])(\.\d+)?$"
)


def looks_like_chembl_id(identifier: str) -> bool:
    return bool(_CHEMBL_ID_RE.match(identifier.strip()))


def looks_like_uniprot_accession(identifier: str) -> bool:
    return bool(_UNIPROT_ACCESSION_RE.match(identifier.strip()))


def entry_name_to_accession(entry_name: str) -> str:
    """UniProt entry name(例: CDK4_HUMAN)をaccession(例: P11802)に変換する。"""
    logger.info("Resolving UniProt entry name %s -> accession ...", entry_name)
    resp = requests.get(
        UNIPROT_SEARCH_URL,
        params={"query": f"id:{entry_name}", "fields": "accession", "format": "json", "size": 1},
        timeout=30,
    )
    resp.raise_for_status()
    results = resp.json().get("results", [])
    if not results:
        raise ValueError(f"UniProt entry name not found: {entry_name}")
    accession = results[0]["primaryAccession"]
    logger.info("  -> %s", accession)
    return accession


def accession_to_chembl_target_id(accession: str) -> str:
    """UniProt accession(例: P11802)をChEMBL target id(例: CHEMBL331)に変換する。"""
    logger.info("Resolving UniProt accession %s -> ChEMBL target id ...", accession)
    resp = requests.get(
        CHEMBL_TARGET_URL,
        params={"target_components__accession": accession, "format": "json"},
        timeout=30,
    )
    resp.raise_for_status()
    targets = resp.json().get("targets", [])
    if not targets:
        raise ValueError(f"ChEMBL target not found for UniProt accession: {accession}")
    chembl_id = targets[0]["target_chembl_id"]
    logger.info("  -> %s", chembl_id)
    return chembl_id


def resolve_uniprot_accession(identifier: str) -> str:
    """UniProt entry name(例: CDK4_HUMAN)/ accession(例: P11802)のいずれを与えても
    UniProt accessionを返す。
    """
    identifier = identifier.strip()
    return identifier if looks_like_uniprot_accession(identifier) else entry_name_to_accession(identifier)


def resolve_chembl_target_id(identifier: str) -> str:
    """UniProt entry name / accession / ChEMBL target idのいずれを与えても
    ChEMBL target idを返す。
    """
    identifier = identifier.strip()
    if looks_like_chembl_id(identifier):
        return identifier.upper()
    return accession_to_chembl_target_id(resolve_uniprot_accession(identifier))

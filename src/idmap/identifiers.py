"""UniProt entry name / accession と ChEMBL target id の相互マッピング。

現時点では「識別子 -> ChEMBL target id」「識別子 -> UniProt accession」の方向のみ実装している。
逆方向の解決が必要になった時点で、この下に追加していく。
"""

import re

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{identifier}.json"
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
    """UniProt entry name(例: CDK4_HUMAN)をaccession(例: P11802)に変換する。

    エントリ取得エンドポイント(`/uniprotkb/{id}.json`)はentry name/accessionのどちらを
    与えても該当accessionへ直接解決される(あいまいさのない一意な解決)。検索エンドポイント
    (`query=id:...`)はトークン化された全文検索であり、無関係のエントリ(例: `CDK1_HUMAN`に対する
    `CDKA1_HUMAN`)が先頭に返ることがあるため使わない。
    """
    entry_name = entry_name.strip()
    logger.info("Resolving UniProt entry name %s -> accession ...", entry_name)
    resp = requests.get(UNIPROT_ENTRY_URL.format(identifier=entry_name), timeout=30)
    if resp.status_code == 400:
        raise ValueError(f"UniProt entry name not found: {entry_name}")
    resp.raise_for_status()
    accession = resp.json()["primaryAccession"]
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

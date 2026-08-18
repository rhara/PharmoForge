"""ローカルのChEMBL SQLiteデータベース(ChEMBL公式配布のchembl_XX.db)から、
ChEMBL Web APIと同等の情報を取得する。

Web API(www.ebi.ac.uk/chembl/api/data/...)が障害・レート制限等で使えない場合のフォールバック。
"""

import sqlite3
from pathlib import Path

from core.logging_utils import get_logger

logger = get_logger(__name__)


def resolve_target_chembl_id(accession: str, db_path: str | Path) -> str:
    """UniProt accessionから、SINGLE PROTEINターゲットのChEMBL target idを解決する。

    `idmap.accession_to_chembl_target_id`(Web API版)と同じ入出力。見つからなければ`ValueError`。
    """
    logger.info("Resolving UniProt accession %s -> ChEMBL target id (local DB) ...", accession)
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.execute(
            """
            SELECT DISTINCT td.chembl_id
            FROM component_sequences cs
            JOIN target_components tc ON cs.component_id = tc.component_id
            JOIN target_dictionary td ON tc.tid = td.tid
            WHERE cs.accession = ? AND td.target_type = 'SINGLE PROTEIN'
            """,
            (accession,),
        )
        rows = cur.fetchall()
    finally:
        conn.close()

    if not rows:
        raise ValueError(f"ChEMBL target not found for UniProt accession: {accession}")
    chembl_id = rows[0][0]
    logger.info("  -> %s", chembl_id)
    return chembl_id


def fetch_activities(target_chembl_id: str, db_path: str | Path) -> list[dict]:
    """指定したChEMBL target idについて、pChEMBL値を持つ活性データを取得する(ローカルDB版)。

    `chembl.activity.fetch_activities`(Web API版)と同じ主要フィールドを持つdictのリストを返す:
    `molecule_chembl_id`/`molecule_pref_name`/`canonical_smiles`/`standard_type`/`standard_value`/
    `standard_units`/`pchembl_value`/`assay_chembl_id`/`document_chembl_id`。
    """
    logger.info(
        "Fetching activities for target %s (pChEMBL value required, local DB) ...", target_chembl_id
    )
    conn = sqlite3.connect(str(db_path))
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.execute(
            """
            SELECT
                md.chembl_id AS molecule_chembl_id,
                md.pref_name AS molecule_pref_name,
                cs.canonical_smiles,
                act.standard_type,
                act.standard_value,
                act.standard_units,
                act.pchembl_value,
                a.chembl_id AS assay_chembl_id,
                d.chembl_id AS document_chembl_id
            FROM activities act
            JOIN assays a ON act.assay_id = a.assay_id
            JOIN target_dictionary td ON a.tid = td.tid
            JOIN molecule_dictionary md ON act.molregno = md.molregno
            LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
            LEFT JOIN docs d ON act.doc_id = d.doc_id
            WHERE td.chembl_id = ? AND act.pchembl_value IS NOT NULL
            """,
            (target_chembl_id,),
        )
        records = [dict(row) for row in cur.fetchall()]
    finally:
        conn.close()

    logger.info("  -> %d activities fetched for %s", len(records), target_chembl_id)
    return records

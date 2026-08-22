"""NCBI BLAST Web API(QBLAST)を用いた配列相同性検索(requestsベース)。"""

import re
import time

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def submit_blast(
    sequence: str,
    program: str = "blastp",
    database: str = "pdb",
    entrez_query: str | None = None,
) -> str:
    """BLASTジョブを投函し、Request ID(RID)を返す。

    `entrez_query`はNCBI Entrezのクエリ構文で検索対象を絞り込む(例: `"Homo sapiens[Organism]"`)。
    """
    logger.info(
        "Submitting %s search against %s database (query length=%d, entrez_query=%s) ...",
        program,
        database,
        len(sequence),
        entrez_query,
    )
    data = {"CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": sequence}
    if entrez_query:
        data["ENTREZ_QUERY"] = entrez_query
    resp = requests.post(NCBI_BLAST_URL, data=data, timeout=60)
    resp.raise_for_status()
    match = re.search(r"RID = (\S+)", resp.text)
    if not match:
        raise RuntimeError("Failed to submit BLAST job: RID not found in response")
    rid = match.group(1)
    logger.info("  -> RID=%s", rid)
    return rid


def wait_for_blast(rid: str, poll_interval: float = 10.0, timeout: float = 600.0) -> None:
    """BLASTジョブの完了(Status=READY)までポーリングして待機する。"""
    logger.info("Waiting for BLAST job %s to complete ...", rid)
    elapsed = 0.0
    while True:
        resp = requests.get(
            NCBI_BLAST_URL,
            params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid},
            timeout=30,
        )
        resp.raise_for_status()
        status_match = re.search(r"Status=(\S+)", resp.text)
        status = status_match.group(1) if status_match else "UNKNOWN"
        logger.info("  [%.0fs] status=%s", elapsed, status)

        if status == "READY":
            if "ThereAreHits=yes" not in resp.text:
                logger.warning("BLAST job %s finished with no hits", rid)
            return
        if status != "WAITING":
            raise RuntimeError(f"BLAST job {rid} failed with status {status}")
        if elapsed >= timeout:
            raise TimeoutError(f"BLAST job {rid} did not complete within {timeout}s")

        time.sleep(poll_interval)
        elapsed += poll_interval


def parse_pdb_subject_id(subject_id: str) -> tuple[str, str]:
    """database="pdb"のヒットのsubject id(例: pdb|6GZM|A, 6GZM_A)からPDB IDとchainを取り出す。"""
    subject_id = subject_id.strip()
    if subject_id.startswith("pdb|"):
        parts = subject_id.split("|")
        return parts[1].upper(), parts[2] if len(parts) > 2 else ""
    if "_" in subject_id:
        pdb_id, chain = subject_id.split("_", 1)
        return pdb_id.upper(), chain
    return subject_id.upper(), ""


def parse_uniprot_subject_id(subject_id: str) -> str:
    """database="swissprot"等のヒットのsubject id(例: sp|P11802|CDK4_HUMAN, P11802.1)から
    UniProt accessionを取り出す(バージョン番号は除く)。
    """
    subject_id = subject_id.strip()
    if "|" in subject_id:
        parts = subject_id.split("|")
        if len(parts) >= 2:
            subject_id = parts[1]
    return subject_id.split(".")[0].upper()


def fetch_hits(rid: str) -> list[dict]:
    """完了したBLASTジョブのヒット一覧をタブ区切り形式で取得しparseする。"""
    logger.info("Fetching BLAST hits for RID=%s ...", rid)
    resp = requests.get(
        NCBI_BLAST_URL,
        params={
            "CMD": "Get",
            "RESULTS_FILE": "on",
            "FORMAT_TYPE": "Text",
            "FORMAT_OBJECT": "Alignment",
            "ALIGNMENT_VIEW": "Tabular",
            "DESCRIPTIONS": 250,
            "RID": rid,
        },
        timeout=60,
    )
    resp.raise_for_status()

    hits: list[dict] = []
    for line in resp.text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        (
            _query_id,
            subject_id,
            pct_identity,
            align_length,
            _mismatches,
            _gap_opens,
            _q_start,
            _q_end,
            _s_start,
            _s_end,
            evalue,
            bit_score,
        ) = fields[:12]
        hits.append(
            {
                "subject_id": subject_id,
                "identity": float(pct_identity),
                "align_length": int(align_length),
                "evalue": float(evalue),
                "bit_score": float(bit_score),
            }
        )
    logger.info("  -> %d hits", len(hits))
    return hits


def format_evalue(evalue: float) -> str:
    """e-valueを表示用に整形する。0はそのまま「0」、それ以外は有効数字2桁の指数表記
    (例: 9.47e-95 -> 9.5e-95)にする。
    """
    return "0" if evalue == 0 else f"{evalue:.1e}"


def best_hit_per_accession(hits: list[dict], exclude_accession: str | None = None) -> list[dict]:
    """UniProt accessionごとに最良ヒット(evalue最小)だけを残し、evalue昇順で返す。

    各ヒットの`subject_id`を`parse_uniprot_subject_id`でaccessionに変換し(`database="swissprot"`等の
    ヒットを想定)、ヒットdictに`"accession"`キーを追加する。`exclude_accession`(クエリ自身の
    accession等)を指定した場合はそのaccessionのヒットを除く。
    """
    best_by_accession: dict[str, dict] = {}
    for h in hits:
        acc = parse_uniprot_subject_id(h["subject_id"])
        if exclude_accession is not None and acc == exclude_accession:
            continue
        if acc not in best_by_accession or h["evalue"] < best_by_accession[acc]["evalue"]:
            best_by_accession[acc] = {**h, "accession": acc}
    return sorted(best_by_accession.values(), key=lambda h: h["evalue"])


def blast_search(
    sequence: str,
    program: str = "blastp",
    database: str = "pdb",
    entrez_query: str | None = None,
    poll_interval: float = 10.0,
    timeout: float = 600.0,
) -> list[dict]:
    """配列を投函し、完了を待ってヒット一覧を返す(submit + wait + fetchの組み合わせ)。"""
    rid = submit_blast(sequence, program=program, database=database, entrez_query=entrez_query)
    wait_for_blast(rid, poll_interval=poll_interval, timeout=timeout)
    return fetch_hits(rid)

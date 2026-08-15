"""AlphaFold DBから蛋白の予測構造ファイルを取得する。"""

from pathlib import Path

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


def fetch_structure(accession: str, output: Path, fmt: str | None = None) -> None:
    """UniProt accession(例: P61626)のAlphaFold DB予測構造をダウンロードしoutputへ保存する。

    fmt(cif/pdb)を省略した場合は、outputの拡張子(.cif / .pdb)から決める。
    どちらもなければcif。
    ダウンロードURLはAlphaFold DBのpredictionエンドポイントから都度解決する
    (バージョン番号をURLに固定で埋め込まない)。
    """
    accession = accession.strip().upper()
    fmt = (fmt or output.suffix.lstrip(".") or "cif").lower()

    logger.info("Resolving AlphaFold DB entry for %s ...", accession)
    api_url = AFDB_API_URL.format(accession=accession)
    resp = requests.get(api_url, timeout=60)
    resp.raise_for_status()
    entries = resp.json()
    if not entries:
        raise ValueError(f"AlphaFold DBに予測構造が見つかりません: {accession}")
    entry = entries[0]

    url_key = "cifUrl" if fmt == "cif" else "pdbUrl"
    structure_url = entry.get(url_key)
    if not structure_url:
        raise ValueError(f"AlphaFold DBのレスポンスに{url_key}がありません: {accession}")

    logger.info(
        "Downloading structure %s (%s, v%s) from AlphaFold DB ...",
        accession,
        fmt,
        entry.get("latestVersion"),
    )
    struct_resp = requests.get(structure_url, timeout=60)
    struct_resp.raise_for_status()

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(struct_resp.content)
    logger.info("Done: saved %s (%d bytes) to %s", accession, len(struct_resp.content), output)


def fetch_structures(accessions: list[str], output_dir: Path, fmt: str) -> None:
    """複数のUniProt accessionをまとめてダウンロードし、`output_dir/<ACCESSION>.<fmt>`として保存する。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    total = len(accessions)
    logger.info("Downloading %d structures (%s) into %s ...", total, fmt, output_dir)
    for i, accession in enumerate(accessions, start=1):
        accession = accession.strip().upper()
        logger.info("[%d/%d] %s", i, total, accession)
        fetch_structure(accession, output_dir / f"{accession}.{fmt}", fmt=fmt)
    logger.info("Done: %d structures saved to %s", total, output_dir)

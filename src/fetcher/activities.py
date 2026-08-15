"""ChEMBL活性データの構造標準化・集約・TSV書き出し(fetcher固有の集計ロジック)。

生の活性データ取得は共通パッケージ[`src/chembl`](../chembl)の`fetch_activities`を使う。
"""

import csv
import statistics
from collections import defaultdict
from pathlib import Path

from molstd import standardize_smiles

from core.logging_utils import get_logger

logger = get_logger(__name__)

AGGREGATED_FIELDS = [
    "smiles",
    "_median",
    "_mean",
    "_sd",
    "_n",
]


def standardize_and_aggregate(records: list[dict]) -> list[dict]:
    """化合物構造を`molstd.standardize_smiles`で標準化し、
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

"""スキャフォールド解析結果のレポート出力(TSV・グリッド画像)。"""

from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from core.logging_utils import get_logger

logger = get_logger(__name__)


def write_summary_tsv(summary: pd.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, sep="\t", index=False)
    logger.info("Wrote scaffold summary (%d scaffolds) to %s", len(summary), output)


def render_scaffold_grid(
    summary: pd.DataFrame,
    output: Path,
    top_n: int = 20,
    ascending: bool = False,
) -> None:
    """enrichmentの上位(ascending=False)または下位(ascending=True)top_n件のスキャフォールドを画像化する。"""
    rows = summary.sort_values("enrichment", ascending=ascending).head(top_n)

    mols = []
    legends = []
    for _, row in rows.iterrows():
        mol = Chem.MolFromSmiles(row["scaffold"])
        if mol is None:
            continue
        mols.append(mol)
        legends.append(
            f"n={int(row['n_total'])} high={int(row['n_high'])} low={int(row['n_low'])}\n"
            f"enrich={row['enrichment']:.2f} med={row['median_activity']:.2f}"
        )

    if not mols:
        logger.warning("No scaffolds to render, skipping %s", output)
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(250, 220), legends=legends, returnPNG=False)
    img.save(output)
    logger.info("Wrote scaffold grid image (%d scaffolds) to %s", len(mols), output)

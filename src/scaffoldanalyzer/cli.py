from pathlib import Path

import click
import pandas as pd

from . import clustering, report

DEFAULT_HIGH_QUANTILE = 0.75
DEFAULT_LOW_QUANTILE = 0.25
DEFAULT_MIN_COUNT = 2
DEFAULT_TOP_N = 20


@click.command("analyze-scaffolds")
@click.argument("input_path", type=click.Path(exists=True, path_type=Path))
@click.option("--smiles-col", default="smiles", show_default=True, help="SMILES column name")
@click.option("--activity-col", default="_median", show_default=True, help="Activity value column name")
@click.option(
    "--output-dir", "-o", required=True, type=click.Path(path_type=Path), help="Output directory"
)
@click.option(
    "--high-quantile",
    default=DEFAULT_HIGH_QUANTILE,
    show_default=True,
    help="Quantile threshold (and above) considered high activity",
)
@click.option(
    "--low-quantile",
    default=DEFAULT_LOW_QUANTILE,
    show_default=True,
    help="Quantile threshold (and below) considered low activity",
)
@click.option(
    "--min-count",
    default=DEFAULT_MIN_COUNT,
    show_default=True,
    help="Minimum occurrence count (summed over all bins) for a scaffold to be included",
)
@click.option(
    "--top-n",
    default=DEFAULT_TOP_N,
    show_default=True,
    help="Number of top/bottom entries drawn in the grid image",
)
def analyze_scaffolds_cmd(
    input_path: Path,
    smiles_col: str,
    activity_col: str,
    output_dir: Path,
    high_quantile: float,
    low_quantile: float,
    min_count: int,
    top_n: int,
):
    """Compare high/low activity distribution by Bemis-Murcko scaffold.

    INPUT_PATH is a TSV/CSV with a SMILES column (default `smiles`) and an activity value
    column (default `_median`); the output of `pf fetch activities=...` can be used as-is.

    \b
    Examples:
      pf analyze-scaffolds data/cdk4_human_activities.tsv --output-dir data/cdk4_scaffold_analysis
    """
    df = pd.read_csv(input_path, sep=None, engine="python")

    df = clustering.add_scaffolds(df, smiles_col=smiles_col)
    df = clustering.assign_activity_bins(
        df, activity_col=activity_col, high_quantile=high_quantile, low_quantile=low_quantile
    )
    summary = clustering.summarize_scaffolds(df, activity_col=activity_col, min_count=min_count)

    output_dir.mkdir(parents=True, exist_ok=True)
    report.write_summary_tsv(summary, output_dir / "scaffold_summary.tsv")
    report.render_scaffold_grid(
        summary, output_dir / "scaffold_grid_high.png", top_n=top_n, ascending=False
    )
    report.render_scaffold_grid(
        summary, output_dir / "scaffold_grid_low.png", top_n=top_n, ascending=True
    )
    report.render_compound_table(
        df,
        summary,
        activity_col,
        output_dir / "scaffold_compounds_high.html",
        top_n=top_n,
        ascending=False,
        smiles_col=smiles_col,
    )
    report.render_compound_table(
        df,
        summary,
        activity_col,
        output_dir / "scaffold_compounds_low.html",
        top_n=top_n,
        ascending=True,
        smiles_col=smiles_col,
    )

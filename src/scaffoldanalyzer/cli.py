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
@click.option("--smiles-col", default="smiles", show_default=True, help="SMILES列名")
@click.option("--activity-col", default="_median", show_default=True, help="活性値列名")
@click.option(
    "--output-dir", "-o", required=True, type=click.Path(path_type=Path), help="出力ディレクトリ"
)
@click.option(
    "--high-quantile", default=DEFAULT_HIGH_QUANTILE, show_default=True, help="高活性とみなす分位点(以上)"
)
@click.option(
    "--low-quantile", default=DEFAULT_LOW_QUANTILE, show_default=True, help="低活性とみなす分位点(以下)"
)
@click.option(
    "--min-count",
    default=DEFAULT_MIN_COUNT,
    show_default=True,
    help="集計対象とするスキャフォールドの最小出現数(全bin合計)",
)
@click.option(
    "--top-n", default=DEFAULT_TOP_N, show_default=True, help="グリッド画像に描画する上位・下位件数"
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
    """Bemis-Murckoスキャフォールドで活性の高低分布を比較する。

    INPUT_PATHはSMILES列(既定`smiles`)と活性値列(既定`_median`)を持つTSV/CSV
    (`pf fetch activities=...`の出力をそのまま使える)。

    \b
    例:
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

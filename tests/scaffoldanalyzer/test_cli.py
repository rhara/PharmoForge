import pandas as pd
from click.testing import CliRunner

from scaffoldanalyzer.cli import analyze_scaffolds_cmd


def test_analyze_scaffolds_cmd_end_to_end(tmp_path):
    input_tsv = tmp_path / "activities.tsv"
    pd.DataFrame(
        {
            "smiles": [
                "Cc1ccccc1",
                "CCc1ccccc1",
                "CCCc1ccccc1",
                "CCO",
                "CCN",
                "CCC",
                "CCCC",
                "CCCCC",
            ],
            "_median": [9.0, 8.8, 8.5, 8.0, 3.0, 2.8, 2.5, 2.0],
        }
    ).to_csv(input_tsv, sep="\t", index=False)

    output_dir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(
        analyze_scaffolds_cmd,
        [str(input_tsv), "--output-dir", str(output_dir), "--min-count", "1"],
    )

    assert result.exit_code == 0, result.output
    assert (output_dir / "scaffold_summary.tsv").exists()
    assert (output_dir / "scaffold_grid_high.png").exists()
    assert (output_dir / "scaffold_grid_low.png").exists()

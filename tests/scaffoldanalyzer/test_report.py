import pandas as pd

from scaffoldanalyzer import report


def _summary_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "scaffold": ["c1ccccc1", "c1ccncc1"],
            "n_total": [3, 2],
            "n_high": [2, 0],
            "n_mid": [1, 1],
            "n_low": [0, 1],
            "frac_high": [0.5, 0.0],
            "frac_low": [0.0, 0.3],
            "enrichment": [0.5, -0.3],
            "mean_activity": [8.5, 3.0],
            "median_activity": [8.5, 3.0],
        }
    )


def test_write_summary_tsv(tmp_path):
    output = tmp_path / "summary.tsv"
    report.write_summary_tsv(_summary_df(), output)

    lines = output.read_text().splitlines()
    assert lines[0].split("\t")[0] == "scaffold"
    assert len(lines) == 3


def test_render_scaffold_grid_creates_image(tmp_path):
    output = tmp_path / "grid.png"
    report.render_scaffold_grid(_summary_df(), output, top_n=5, ascending=False)

    assert output.exists()
    assert output.stat().st_size > 0


def test_render_scaffold_grid_skips_when_empty(tmp_path):
    output = tmp_path / "grid.png"
    empty = _summary_df().iloc[0:0]
    report.render_scaffold_grid(empty, output, top_n=5)

    assert not output.exists()

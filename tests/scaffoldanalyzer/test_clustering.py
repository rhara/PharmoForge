import pandas as pd

from scaffoldanalyzer.clustering import add_scaffolds, assign_activity_bins, summarize_scaffolds


def _sample_df() -> pd.DataFrame:
    # トルエン誘導体3件(高活性寄り) + 無関係な低活性化合物2件
    return pd.DataFrame(
        {
            "smiles": ["Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1", "CCO", "not a smiles"],
            "activity": [9.0, 8.5, 8.0, 2.0, 1.0],
        }
    )


def test_add_scaffolds_drops_unparsable_rows():
    df = add_scaffolds(_sample_df(), smiles_col="smiles")
    assert len(df) == 4
    assert (df["scaffold"] == "c1ccccc1").sum() == 3


def test_summarize_scaffolds_enrichment():
    df = add_scaffolds(_sample_df(), smiles_col="smiles")
    df = assign_activity_bins(df, "activity", high_quantile=0.6, low_quantile=0.4)
    summary = summarize_scaffolds(df, activity_col="activity", min_count=1)

    benzene_row = summary[summary["scaffold"] == "c1ccccc1"].iloc[0]
    assert benzene_row["n_total"] == 3
    assert benzene_row["n_high"] >= 1
    assert benzene_row["enrichment"] > 0
    # enrichment降順でソートされている
    assert summary["enrichment"].is_monotonic_decreasing


def test_summarize_scaffolds_columns_and_rounding():
    df = add_scaffolds(_sample_df(), smiles_col="smiles")
    df = assign_activity_bins(df, "activity", high_quantile=0.6, low_quantile=0.4)
    summary = summarize_scaffolds(df, activity_col="activity", min_count=1)

    assert list(summary.columns) == [
        "scaffold",
        "n_total",
        "n_high",
        "n_mid",
        "n_low",
        "frac_high",
        "frac_low",
        "enrichment",
        "median_activity",
        "mean_activity",
        "sd_activity",
    ]

    benzene_row = summary[summary["scaffold"] == "c1ccccc1"].iloc[0]
    assert benzene_row["median_activity"] == 8.5
    assert benzene_row["mean_activity"] == 8.5
    assert round(benzene_row["sd_activity"], 2) == benzene_row["sd_activity"]

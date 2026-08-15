from unittest.mock import patch

from fetcher import activities


@patch("fetcher.activities.standardize_smiles")
def test_standardize_and_aggregate_groups_by_standardized_structure(mock_standardize):
    # 2 records standardize to the same structure ("CCO"), 1 to a different one ("CCN")
    mock_standardize.side_effect = lambda smiles: {"CCO.raw": "CCO", "CCO.raw2": "CCO", "CCN.raw": "CCN"}[smiles]
    records = [
        {"canonical_smiles": "CCO.raw", "pchembl_value": "6.0"},
        {"canonical_smiles": "CCO.raw2", "pchembl_value": "8.0"},
        {"canonical_smiles": "CCN.raw", "pchembl_value": "5.0"},
    ]

    aggregated = activities.standardize_and_aggregate(records)
    by_smiles = {row["smiles"]: row for row in aggregated}

    assert set(by_smiles) == {"CCO", "CCN"}
    assert by_smiles["CCO"]["_n"] == 2
    assert by_smiles["CCO"]["_mean"] == 7.0
    assert by_smiles["CCO"]["_median"] == 7.0
    assert by_smiles["CCO"]["_sd"] == round(((6.0 - 7.0) ** 2 + (8.0 - 7.0) ** 2) ** 0.5, 3)
    assert by_smiles["CCN"]["_n"] == 1
    assert by_smiles["CCN"]["_sd"] == ""


@patch("fetcher.activities.standardize_smiles")
def test_standardize_and_aggregate_sorted_by_median_descending(mock_standardize):
    mock_standardize.side_effect = lambda smiles: smiles
    records = [
        {"canonical_smiles": "LOW", "pchembl_value": "5.0"},
        {"canonical_smiles": "HIGH", "pchembl_value": "9.0"},
        {"canonical_smiles": "MID", "pchembl_value": "7.0"},
    ]

    aggregated = activities.standardize_and_aggregate(records)

    assert [row["smiles"] for row in aggregated] == ["HIGH", "MID", "LOW"]


@patch("fetcher.activities.standardize_smiles", return_value=None)
def test_standardize_and_aggregate_skips_unparsable(mock_standardize):
    records = [{"canonical_smiles": "garbage", "pchembl_value": "6.0"}]
    assert activities.standardize_and_aggregate(records) == []


def test_standardize_and_aggregate_skips_missing_fields():
    records = [
        {"canonical_smiles": None, "pchembl_value": "6.0"},
        {"canonical_smiles": "CCO", "pchembl_value": None},
    ]
    assert activities.standardize_and_aggregate(records) == []


def test_write_activities_tsv(tmp_path):
    records = [
        {
            "smiles": "CCO",
            "_mean": 7.0,
            "_median": 7.0,
            "_sd": 1.414,
            "_n": 2,
        }
    ]
    output = tmp_path / "out.tsv"

    activities.write_activities_tsv(records, output)

    lines = output.read_text().splitlines()
    assert lines[0] == "\t".join(activities.AGGREGATED_FIELDS)
    assert "CCO" in lines[1]
    assert "7.0" in lines[1]

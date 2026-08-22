from unittest.mock import patch

import numpy as np
import pandas as pd

from chembl.aggregate import (
    collect_standardized_activities,
    rollup_compound_summary,
    select_high_potency_compounds,
    summarize_compound_target_activity,
)


def test_collect_standardized_activities_skips_missing_target_and_bad_records():
    targets_df = pd.DataFrame([
        {"accession": "P24941", "entry_name": "CDK2_HUMAN", "chembl_target_id": "CHEMBL301"},
        {"accession": "P99999", "entry_name": "NOTARGET_HUMAN", "chembl_target_id": np.nan},
    ])
    fake_activities = [
        {"canonical_smiles": "c1ccccc1", "pchembl_value": 8.0},
        {"canonical_smiles": None, "pchembl_value": 7.0},  # SMILES欠損は除外
        {"canonical_smiles": "c1ccccc1", "pchembl_value": None},  # pChEMBL欠損は除外
        {"canonical_smiles": "bad smiles", "pchembl_value": 6.0},  # 標準化失敗は除外
    ]

    def fake_standardize(smiles):
        return "c1ccccc1" if smiles == "c1ccccc1" else None

    with (
        patch("chembl.aggregate.fetch_activities", return_value=fake_activities) as mock_fetch,
        patch("chembl.aggregate.standardize_smiles", side_effect=fake_standardize),
    ):
        result = collect_standardized_activities(targets_df, "dummy.db")

    mock_fetch.assert_called_once_with("CHEMBL301", "dummy.db")  # NaN target_idの行はfetchされない
    assert len(result) == 1
    assert result.iloc[0]["smiles"] == "c1ccccc1"
    assert result.iloc[0]["accession"] == "P24941"
    assert result.iloc[0]["pchembl_value"] == 8.0


def test_summarize_compound_target_activity_aggregates_median_sorted_desc():
    activity_df = pd.DataFrame([
        {"smiles": "A", "accession": "P1", "entry_name": "T1", "pchembl_value": 7.0},
        {"smiles": "A", "accession": "P1", "entry_name": "T1", "pchembl_value": 9.0},
        {"smiles": "B", "accession": "P1", "entry_name": "T1", "pchembl_value": 5.0},
    ])

    result = summarize_compound_target_activity(activity_df)

    assert list(result["smiles"]) == ["A", "B"]  # median降順(Aのmedian=8.0 > Bのmedian=5.0)
    assert result.iloc[0]["median"] == 8.0
    assert result.iloc[0]["count"] == 2


def test_rollup_compound_summary_picks_best_target_and_counts():
    activity_summary_df = pd.DataFrame([
        {"smiles": "A", "accession": "P1", "entry_name": "T1", "median": 7.0},
        {"smiles": "A", "accession": "P2", "entry_name": "T2", "median": 9.0},
        {"smiles": "B", "accession": "P1", "entry_name": "T1", "median": 6.0},
    ])

    result = rollup_compound_summary(activity_summary_df)

    a_row = result[result["smiles"] == "A"].iloc[0]
    assert a_row["best_target_entry_name"] == "T2"
    assert a_row["best_median"] == 9.0
    assert a_row["target_count"] == 2

    b_row = result[result["smiles"] == "B"].iloc[0]
    assert b_row["target_count"] == 1
    assert list(result["smiles"]) == ["A", "B"]  # best_median降順


def test_select_high_potency_compounds_filters_potency_only():
    df = pd.DataFrame([
        {"smiles": "A", "best_pchembl_median": 9.5},
        {"smiles": "B", "best_pchembl_median": 8.0},
    ])

    result = select_high_potency_compounds(df, potency_col="best_pchembl_median", potency_cutoff=9.0)

    assert list(result["smiles"]) == ["A"]


def test_select_high_potency_compounds_filters_potency_and_mol_weight():
    df = pd.DataFrame([
        {"smiles": "A", "best_pchembl_median": 9.5, "mol_weight": 300.0},
        {"smiles": "B", "best_pchembl_median": 9.5, "mol_weight": 900.0},  # 分子量オーバー
        {"smiles": "C", "best_pchembl_median": 8.0, "mol_weight": 300.0},  # 活性不足
    ])

    result = select_high_potency_compounds(
        df, potency_col="best_pchembl_median", potency_cutoff=9.0,
        mol_weight_col="mol_weight", mol_weight_range=(250, 650),
    )

    assert list(result["smiles"]) == ["A"]

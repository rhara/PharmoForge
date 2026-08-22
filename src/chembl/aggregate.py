"""複数のChEMBL標的にまたがる活性データの集計(化合物×標的単位の要約、化合物単位のロールアップ)。"""

from pathlib import Path

import pandas as pd
from molstd import standardize_smiles

from core.logging_utils import get_logger

from .local import fetch_activities

logger = get_logger(__name__)


def collect_standardized_activities(
    targets_df: pd.DataFrame,
    db_path: str | Path,
    target_id_col: str = "chembl_target_id",
    accession_col: str = "accession",
    entry_name_col: str = "entry_name",
) -> pd.DataFrame:
    """複数のChEMBL標的について活性データを収集し、化合物ごとの標準化SMILES・pChEMBL値のみの
    フラットなテーブル(列: `smiles`/`accession`/`entry_name`/`pchembl_value`)にする
    (ローカルDB版、`chembl.local.fetch_activities`を使う)。

    `target_id_col`がNaN(ChEMBL targetが見つからなかった標的)の行はスキップする。pChEMBL値
    またはSMILESを欠く活性、標準化に失敗した化合物(`molstd.standardize_smiles`がNoneを返すもの)
    は除外する。
    """
    targets = targets_df.dropna(subset=[target_id_col])
    records = []
    total = len(targets)
    for i, row in enumerate(targets.itertuples(), start=1):
        target_id = getattr(row, target_id_col)
        accession = getattr(row, accession_col)
        entry_name = getattr(row, entry_name_col)
        activities = fetch_activities(target_id, db_path)
        logger.info(
            "[%d/%d] %s (%s): %d activities, standardizing ...", i, total, entry_name, target_id, len(activities)
        )
        n_ok = 0
        for a in activities:
            if a["pchembl_value"] is None or a["canonical_smiles"] is None:
                continue
            std_smiles = standardize_smiles(a["canonical_smiles"])
            if std_smiles is None:
                continue
            records.append({
                "smiles": std_smiles,
                "accession": accession,
                "entry_name": entry_name,
                "pchembl_value": float(a["pchembl_value"]),
            })
            n_ok += 1
        logger.info("    -> %d usable records", n_ok)
    return pd.DataFrame(records, columns=["smiles", "accession", "entry_name", "pchembl_value"])


def summarize_compound_target_activity(
    activity_df: pd.DataFrame,
    group_cols: tuple[str, ...] = ("smiles", "accession", "entry_name"),
    value_col: str = "pchembl_value",
) -> pd.DataFrame:
    """化合物×標的の単位で活性値のmedian/mean/std/個数を集計する(median降順)。

    同じ化合物(標準化SMILES)が同じ標的に対して複数の活性データを持つ場合(異なるアッセイ等)に、
    median(代表値、外れ値に頑健)を主指標として使う想定。
    """
    return (
        activity_df
        .groupby(list(group_cols))[value_col]
        .agg(median="median", mean="mean", std="std", count="count")
        .reset_index()
        .sort_values("median", ascending=False)
        .reset_index(drop=True)
    )


def rollup_compound_summary(
    activity_summary_df: pd.DataFrame,
    smiles_col: str = "smiles",
    entry_name_col: str = "entry_name",
    value_col: str = "median",
) -> pd.DataFrame:
    """`summarize_compound_target_activity`の出力(化合物×標的単位)を、化合物単位にロールアップする。

    同じ化合物が複数の標的に対して活性データを持つ場合があるため、化合物ごとにグループ化し、
    テストされた標的数(`target_count`)・最も高い活性を示した標的名(`best_target_<entry_name_col>`)・
    その活性値(`best_<value_col>`)にまとめる(`best_<value_col>`降順)。
    """
    target_counts = activity_summary_df.groupby(smiles_col).size().rename("target_count").reset_index()

    best_col = f"best_{value_col}"
    best_target_col = f"best_target_{entry_name_col}"
    best_rows = activity_summary_df.loc[activity_summary_df.groupby(smiles_col)[value_col].idxmax()][
        [smiles_col, entry_name_col, value_col]
    ].rename(columns={entry_name_col: best_target_col, value_col: best_col})

    return best_rows.merge(target_counts, on=smiles_col).sort_values(best_col, ascending=False).reset_index(drop=True)


def select_high_potency_compounds(
    df: pd.DataFrame,
    potency_col: str,
    potency_cutoff: float,
    mol_weight_col: str | None = None,
    mol_weight_range: tuple[float, float] | None = None,
) -> pd.DataFrame:
    """`potency_col >= potency_cutoff`の行を抽出する。

    `mol_weight_col`と`mol_weight_range`(`(min, max)`)を両方指定した場合は、分子量でも絞り込む
    (創薬で一般的な「drug-like」範囲によるフィルタ、例: 250〜650)。
    """
    mask = df[potency_col] >= potency_cutoff
    if mol_weight_col is not None and mol_weight_range is not None:
        low, high = mol_weight_range
        mask &= df[mol_weight_col].between(low, high)
    result = df[mask].reset_index(drop=True)
    logger.info("%d/%d compounds pass the potency/mol_weight filter", len(result), len(df))
    return result

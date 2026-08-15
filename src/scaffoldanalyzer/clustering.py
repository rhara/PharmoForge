"""Bemis-Murckoスキャフォールドに基づく活性データのクラスタリング解析。"""

import pandas as pd
from actbin import assign_activity_bins
from molscaffold import compute_scaffold

from core.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = ["add_scaffolds", "assign_activity_bins", "summarize_scaffolds"]


def add_scaffolds(df: pd.DataFrame, smiles_col: str = "smiles") -> pd.DataFrame:
    """各化合物にスキャフォールド列(`scaffold`)を追加する。パース失敗行は除外する。"""
    scaffolds = df[smiles_col].map(compute_scaffold)
    n_failed = int(scaffolds.isna().sum())
    if n_failed:
        logger.warning("Failed to compute scaffold for %d/%d compounds (excluded)", n_failed, len(df))
    return df.assign(scaffold=scaffolds).dropna(subset=["scaffold"]).reset_index(drop=True)


def summarize_scaffolds(
    df: pd.DataFrame,
    activity_col: str,
    min_count: int = 2,
) -> pd.DataFrame:
    """スキャフォールドごとにhigh/mid/low件数とactivity統計量を集計する。

    `df`は`scaffold`列(add_scaffolds)と`bin`列(assign_activity_bins)を持つこと。
    `enrichment`(high群での出現割合 - low群での出現割合)の降順でソートして返す。
    出現総数が`min_count`未満のスキャフォールドは、解釈の信頼性が低いため除外する。
    """
    counts = df.groupby(["scaffold", "bin"]).size().unstack(fill_value=0)
    for col in ("high", "mid", "low"):
        if col not in counts.columns:
            counts[col] = 0
    counts["n_total"] = counts[["high", "mid", "low"]].sum(axis=1)

    activity_stats = df.groupby("scaffold")[activity_col].agg(["mean", "median"])
    summary = counts.join(activity_stats)
    summary = summary[summary["n_total"] >= min_count]

    n_high_total = int((df["bin"] == "high").sum())
    n_low_total = int((df["bin"] == "low").sum())
    summary["frac_high"] = summary["high"] / n_high_total if n_high_total else 0.0
    summary["frac_low"] = summary["low"] / n_low_total if n_low_total else 0.0
    summary["enrichment"] = summary["frac_high"] - summary["frac_low"]

    summary = summary.rename(
        columns={
            "high": "n_high",
            "mid": "n_mid",
            "low": "n_low",
            "mean": "mean_activity",
            "median": "median_activity",
        }
    )
    summary = summary.reset_index().sort_values("enrichment", ascending=False).reset_index(drop=True)
    return summary[
        [
            "scaffold",
            "n_total",
            "n_high",
            "n_mid",
            "n_low",
            "frac_high",
            "frac_low",
            "enrichment",
            "mean_activity",
            "median_activity",
        ]
    ]

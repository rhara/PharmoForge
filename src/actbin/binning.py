"""活性値(または任意の連続値)を分位点でhigh/mid/lowに分類する。"""

import pandas as pd

from core.logging_utils import get_logger

logger = get_logger(__name__)


def assign_activity_bins(
    df: pd.DataFrame,
    activity_col: str,
    high_quantile: float = 0.75,
    low_quantile: float = 0.25,
) -> pd.DataFrame:
    """activity_colの分位点で各行をhigh/mid/lowに分類する(`bin`列を追加)。

    activity_col >= high_quantile分位点 → high、<= low_quantile分位点 → low、それ以外 → mid。
    """
    if not 0 < low_quantile < high_quantile < 1:
        raise ValueError(
            f"0 < low_quantile < high_quantile < 1 を満たすように指定してください "
            f"(low_quantile={low_quantile}, high_quantile={high_quantile})"
        )
    high_threshold = df[activity_col].quantile(high_quantile)
    low_threshold = df[activity_col].quantile(low_quantile)
    logger.info(
        "Activity bin thresholds (%s): low<=%.3f (q=%.2f), high>=%.3f (q=%.2f)",
        activity_col,
        low_threshold,
        low_quantile,
        high_threshold,
        high_quantile,
    )

    def _bin(value: float) -> str:
        if value >= high_threshold:
            return "high"
        if value <= low_threshold:
            return "low"
        return "mid"

    return df.assign(bin=df[activity_col].map(_bin))

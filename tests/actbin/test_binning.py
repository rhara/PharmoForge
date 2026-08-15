import pandas as pd
import pytest

from actbin import assign_activity_bins


def test_assign_activity_bins_high_mid_low():
    df = pd.DataFrame({"value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]})

    result = assign_activity_bins(df, "value", high_quantile=0.75, low_quantile=0.25)

    assert result.loc[result["value"] <= 3, "bin"].eq("low").all()
    assert result.loc[result["value"] >= 8, "bin"].eq("high").all()
    assert result.loc[(result["value"] > 3) & (result["value"] < 8), "bin"].eq("mid").all()


def test_assign_activity_bins_rejects_invalid_quantile_order():
    df = pd.DataFrame({"value": [1, 2, 3]})
    with pytest.raises(ValueError):
        assign_activity_bins(df, "value", high_quantile=0.25, low_quantile=0.75)

from unittest.mock import MagicMock, patch

from chembl import activity


@patch("chembl.activity.requests.get")
def test_fetch_activities_paginates(mock_get):
    page1 = MagicMock(
        json=lambda: {
            "activities": [{"molecule_chembl_id": "CHEMBL1", "pchembl_value": "6.5"}],
            "page_meta": {"next": "/chembl/api/data/activity.json?offset=1000"},
        },
        raise_for_status=lambda: None,
    )
    page2 = MagicMock(
        json=lambda: {
            "activities": [{"molecule_chembl_id": "CHEMBL2", "pchembl_value": "7.1"}],
            "page_meta": {"next": None},
        },
        raise_for_status=lambda: None,
    )
    mock_get.side_effect = [page1, page2]

    records = activity.fetch_activities("CHEMBL331")

    assert [r["molecule_chembl_id"] for r in records] == ["CHEMBL1", "CHEMBL2"]
    assert mock_get.call_count == 2

from unittest.mock import MagicMock, patch

from rcsb import metadata


@patch("rcsb.metadata.requests.get")
def test_fetch_entry_info_xray_with_resolution(mock_get):
    mock_get.return_value = MagicMock(raise_for_status=lambda: None)
    mock_get.return_value.json.return_value = {
        "rcsb_entry_info": {
            "experimental_method": "X-RAY DIFFRACTION",
            "resolution_combined": [1.8],
        }
    }

    info = metadata.fetch_entry_info("6gzm")

    assert info == {"pdb_id": "6GZM", "method": "X-RAY DIFFRACTION", "resolution": 1.8}
    called_url = mock_get.call_args[0][0]
    assert called_url == "https://data.rcsb.org/rest/v1/core/entry/6GZM"


@patch("rcsb.metadata.requests.get")
def test_fetch_entry_info_no_resolution(mock_get):
    mock_get.return_value = MagicMock(raise_for_status=lambda: None)
    mock_get.return_value.json.return_value = {
        "rcsb_entry_info": {"experimental_method": "SOLUTION NMR"}
    }

    info = metadata.fetch_entry_info("2ABC")

    assert info == {"pdb_id": "2ABC", "method": "SOLUTION NMR", "resolution": None}


@patch("rcsb.metadata.fetch_entry_info")
def test_fetch_entries_info_aggregates(mock_fetch_one):
    mock_fetch_one.side_effect = [
        {"pdb_id": "6GZM", "method": "X-RAY DIFFRACTION", "resolution": 1.8},
        {"pdb_id": "7SJ3", "method": "X-RAY DIFFRACTION", "resolution": 2.1},
    ]

    results = metadata.fetch_entries_info(["6gzm", "7sj3"])

    assert results == [
        {"pdb_id": "6GZM", "method": "X-RAY DIFFRACTION", "resolution": 1.8},
        {"pdb_id": "7SJ3", "method": "X-RAY DIFFRACTION", "resolution": 2.1},
    ]

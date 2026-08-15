from unittest.mock import MagicMock, patch

from fetcher import chembl


@patch("fetcher.chembl.requests.get")
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

    records = chembl.fetch_activities("CHEMBL331")

    assert [r["molecule_chembl_id"] for r in records] == ["CHEMBL1", "CHEMBL2"]
    assert mock_get.call_count == 2


def test_write_activities_tsv(tmp_path):
    records = [
        {
            "molecule_chembl_id": "CHEMBL1",
            "canonical_smiles": "CCO",
            "pchembl_value": "6.5",
        }
    ]
    output = tmp_path / "out.tsv"

    chembl.write_activities_tsv(records, output)

    lines = output.read_text().splitlines()
    assert lines[0] == "\t".join(chembl.ACTIVITY_FIELDS)
    assert "CHEMBL1" in lines[1]
    assert "CCO" in lines[1]

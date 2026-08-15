from unittest.mock import MagicMock, patch

from idmap import identifiers


def test_looks_like_chembl_id():
    assert identifiers.looks_like_chembl_id("CHEMBL331")
    assert not identifiers.looks_like_chembl_id("P11802")


def test_looks_like_uniprot_accession():
    assert identifiers.looks_like_uniprot_accession("P11802")
    assert not identifiers.looks_like_uniprot_accession("CDK4_HUMAN")
    assert not identifiers.looks_like_uniprot_accession("CHEMBL331")


@patch("idmap.identifiers.requests.get")
def test_entry_name_to_accession(mock_get):
    mock_get.return_value = MagicMock(
        json=lambda: {"results": [{"primaryAccession": "P11802"}]},
        raise_for_status=lambda: None,
    )
    assert identifiers.entry_name_to_accession("CDK4_HUMAN") == "P11802"


@patch("idmap.identifiers.requests.get")
def test_accession_to_chembl_target_id(mock_get):
    mock_get.return_value = MagicMock(
        json=lambda: {"targets": [{"target_chembl_id": "CHEMBL331"}]},
        raise_for_status=lambda: None,
    )
    assert identifiers.accession_to_chembl_target_id("P11802") == "CHEMBL331"


def test_resolve_chembl_target_id_passthrough():
    assert identifiers.resolve_chembl_target_id("chembl331") == "CHEMBL331"


@patch("idmap.identifiers.accession_to_chembl_target_id", return_value="CHEMBL331")
def test_resolve_chembl_target_id_from_accession(mock_acc):
    assert identifiers.resolve_chembl_target_id("P11802") == "CHEMBL331"


@patch("idmap.identifiers.accession_to_chembl_target_id", return_value="CHEMBL331")
@patch("idmap.identifiers.entry_name_to_accession", return_value="P11802")
def test_resolve_chembl_target_id_from_entry_name(mock_entry, mock_acc):
    assert identifiers.resolve_chembl_target_id("CDK4_HUMAN") == "CHEMBL331"

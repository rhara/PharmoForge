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
        status_code=200,
        json=lambda: {"primaryAccession": "P11802", "uniProtkbId": "CDK4_HUMAN"},
        raise_for_status=lambda: None,
    )
    assert identifiers.entry_name_to_accession("CDK4_HUMAN") == "P11802"
    called_url = mock_get.call_args[0][0]
    assert called_url == "https://rest.uniprot.org/uniprotkb/CDK4_HUMAN.json"


@patch("idmap.identifiers.requests.get")
def test_entry_name_to_accession_not_found(mock_get):
    mock_get.return_value = MagicMock(status_code=400)
    try:
        identifiers.entry_name_to_accession("NOT_A_REAL_PROTEIN_HUMAN")
        assert False, "expected ValueError"
    except ValueError:
        pass


@patch("idmap.identifiers.requests.get")
def test_accession_to_chembl_target_id(mock_get):
    mock_get.return_value = MagicMock(
        json=lambda: {"targets": [{"target_chembl_id": "CHEMBL331"}]},
        raise_for_status=lambda: None,
    )
    assert identifiers.accession_to_chembl_target_id("P11802") == "CHEMBL331"


def test_resolve_uniprot_accession_passthrough():
    assert identifiers.resolve_uniprot_accession("P11802") == "P11802"


@patch("idmap.identifiers.entry_name_to_accession", return_value="P11802")
def test_resolve_uniprot_accession_from_entry_name(mock_entry):
    assert identifiers.resolve_uniprot_accession("CDK4_HUMAN") == "P11802"
    mock_entry.assert_called_once_with("CDK4_HUMAN")


def test_resolve_chembl_target_id_passthrough():
    assert identifiers.resolve_chembl_target_id("chembl331") == "CHEMBL331"


@patch("idmap.identifiers.accession_to_chembl_target_id", return_value="CHEMBL331")
def test_resolve_chembl_target_id_from_accession(mock_acc):
    assert identifiers.resolve_chembl_target_id("P11802") == "CHEMBL331"


@patch("idmap.identifiers.accession_to_chembl_target_id", return_value="CHEMBL331")
@patch("idmap.identifiers.entry_name_to_accession", return_value="P11802")
def test_resolve_chembl_target_id_from_entry_name(mock_entry, mock_acc):
    assert identifiers.resolve_chembl_target_id("CDK4_HUMAN") == "CHEMBL331"

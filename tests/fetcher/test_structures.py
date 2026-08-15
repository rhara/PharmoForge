from unittest.mock import MagicMock, patch

from fetcher import structures


@patch("fetcher.structures.requests.get")
def test_fetch_structure_cif(mock_get, tmp_path):
    mock_get.return_value = MagicMock(content=b"data_9CSK\n...", raise_for_status=lambda: None)
    output = tmp_path / "9csk.cif"

    structures.fetch_structure("9csk", output)

    assert output.read_bytes() == b"data_9CSK\n..."
    called_url = mock_get.call_args[0][0]
    assert called_url == "https://files.rcsb.org/download/9CSK.cif"


@patch("fetcher.structures.requests.get")
def test_fetch_structure_pdb_format_from_extension(mock_get, tmp_path):
    mock_get.return_value = MagicMock(content=b"HEADER", raise_for_status=lambda: None)
    output = tmp_path / "9csk.pdb"

    structures.fetch_structure("9csk", output)

    called_url = mock_get.call_args[0][0]
    assert called_url == "https://files.rcsb.org/download/9CSK.pdb"

from unittest.mock import MagicMock, patch

from fetcher import structures


@patch("fetcher.structures.requests.get")
def test_fetch_structure_infers_format_from_extension(mock_get, tmp_path):
    mock_get.return_value = MagicMock(content=b"data_9CSK\n...", raise_for_status=lambda: None)
    output = tmp_path / "9csk.cif"

    structures.fetch_structure("9csk", output)

    assert output.read_bytes() == b"data_9CSK\n..."
    called_url = mock_get.call_args[0][0]
    assert called_url == "https://files.rcsb.org/download/9CSK.cif"


@patch("fetcher.structures.requests.get")
def test_fetch_structure_explicit_fmt_overrides_extension(mock_get, tmp_path):
    mock_get.return_value = MagicMock(content=b"HEADER", raise_for_status=lambda: None)
    output = tmp_path / "9csk.cif"

    structures.fetch_structure("9csk", output, fmt="pdb")

    called_url = mock_get.call_args[0][0]
    assert called_url == "https://files.rcsb.org/download/9CSK.pdb"


@patch("fetcher.structures.requests.get")
def test_fetch_structure_pdb_format_from_extension(mock_get, tmp_path):
    mock_get.return_value = MagicMock(content=b"HEADER", raise_for_status=lambda: None)
    output = tmp_path / "9csk.pdb"

    structures.fetch_structure("9csk", output)

    called_url = mock_get.call_args[0][0]
    assert called_url == "https://files.rcsb.org/download/9CSK.pdb"


@patch("fetcher.structures.requests.get")
def test_fetch_structures_downloads_each_into_output_dir(mock_get, tmp_path):
    mock_get.return_value = MagicMock(content=b"data", raise_for_status=lambda: None)
    output_dir = tmp_path / "data"

    structures.fetch_structures(["6p8f", "7sj3", "9csk"], output_dir, "pdb")

    called_urls = [c.args[0] for c in mock_get.call_args_list]
    assert called_urls == [
        "https://files.rcsb.org/download/6P8F.pdb",
        "https://files.rcsb.org/download/7SJ3.pdb",
        "https://files.rcsb.org/download/9CSK.pdb",
    ]
    assert (output_dir / "6P8F.pdb").read_bytes() == b"data"
    assert (output_dir / "7SJ3.pdb").read_bytes() == b"data"
    assert (output_dir / "9CSK.pdb").read_bytes() == b"data"

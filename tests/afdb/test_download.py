from unittest.mock import MagicMock, patch

from afdb import download

API_RESPONSE = [
    {
        "latestVersion": 6,
        "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P61626-F1-model_v6.cif",
        "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P61626-F1-model_v6.pdb",
    }
]


@patch("afdb.download.requests.get")
def test_fetch_structure_infers_format_from_extension(mock_get, tmp_path):
    api_resp = MagicMock(raise_for_status=lambda: None)
    api_resp.json.return_value = API_RESPONSE
    struct_resp = MagicMock(content=b"data_P61626\n...", raise_for_status=lambda: None)
    mock_get.side_effect = [api_resp, struct_resp]
    output = tmp_path / "p61626.cif"

    download.fetch_structure("p61626", output)

    assert output.read_bytes() == b"data_P61626\n..."
    api_url = mock_get.call_args_list[0].args[0]
    struct_url = mock_get.call_args_list[1].args[0]
    assert api_url == "https://alphafold.ebi.ac.uk/api/prediction/P61626"
    assert struct_url == API_RESPONSE[0]["cifUrl"]


@patch("afdb.download.requests.get")
def test_fetch_structure_explicit_fmt_overrides_extension(mock_get, tmp_path):
    api_resp = MagicMock(raise_for_status=lambda: None)
    api_resp.json.return_value = API_RESPONSE
    struct_resp = MagicMock(content=b"HEADER", raise_for_status=lambda: None)
    mock_get.side_effect = [api_resp, struct_resp]
    output = tmp_path / "p61626.cif"

    download.fetch_structure("p61626", output, fmt="pdb")

    struct_url = mock_get.call_args_list[1].args[0]
    assert struct_url == API_RESPONSE[0]["pdbUrl"]


@patch("afdb.download.requests.get")
def test_fetch_structure_raises_when_no_entry(mock_get, tmp_path):
    api_resp = MagicMock(raise_for_status=lambda: None)
    api_resp.json.return_value = []
    mock_get.side_effect = [api_resp]
    output = tmp_path / "unknown.cif"

    try:
        download.fetch_structure("unknown", output)
        assert False, "expected ValueError"
    except ValueError:
        pass


@patch("afdb.download.requests.get")
def test_fetch_structures_downloads_each_into_output_dir(mock_get, tmp_path):
    api_resp = MagicMock(raise_for_status=lambda: None)
    api_resp.json.return_value = API_RESPONSE
    struct_resp = MagicMock(content=b"data", raise_for_status=lambda: None)
    mock_get.side_effect = [api_resp, struct_resp, api_resp, struct_resp]
    output_dir = tmp_path / "data"

    download.fetch_structures(["p61626", "p11802"], output_dir, "pdb")

    assert (output_dir / "P61626.pdb").read_bytes() == b"data"
    assert (output_dir / "P11802.pdb").read_bytes() == b"data"

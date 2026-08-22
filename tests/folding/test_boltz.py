from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from folding.boltz import StructureTemplate, predict_structure, search_msa


def test_resolve_api_key_raises_without_env_or_argument(monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    with pytest.raises(ValueError):
        predict_structure("MKV", "/tmp/out", "x")


def test_search_msa_writes_a3m_file(tmp_path, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "nvapi-test")
    output_path = tmp_path / "seq.a3m"

    with patch("folding.boltz.MSASearchIntegration") as mock_integration_cls:
        mock_integration = mock_integration_cls.return_value
        mock_integration.search_and_save = AsyncMock(return_value=output_path)

        result = search_msa("MKV...", output_path)

    assert result == output_path
    mock_integration.search_and_save.assert_awaited_once()
    call_kwargs = mock_integration.search_and_save.call_args.kwargs
    assert call_kwargs["sequence"] == "MKV..."
    assert call_kwargs["output_format"] == "a3m"


def test_predict_structure_saves_cif_files(tmp_path, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "nvapi-test")
    msa_path = tmp_path / "seq.a3m"
    msa_path.write_text(">Query\nMKV\n")

    fake_response = MagicMock()
    fake_response.structures = [MagicMock(structure="data_model\n#\n")]
    fake_response.confidence_scores = [0.87]
    fake_response.ptm_scores = [0.75]

    with patch("folding.boltz.Boltz2Client") as mock_client_cls:
        mock_client = mock_client_cls.return_value
        mock_client.predict = AsyncMock(return_value=fake_response)

        result = predict_structure(
            "MKV", tmp_path / "out", "test_protein", msa_path=msa_path,
        )

    assert len(result.structure_paths) == 1
    assert result.structure_paths[0].read_text() == "data_model\n#\n"
    assert result.confidence_scores == [0.87]
    assert result.ptm_scores == [0.75]


def test_predict_structure_passes_template_to_polymer(tmp_path, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "nvapi-test")
    template_path = tmp_path / "template.cif"
    template_path.write_text("data_template\n#\n")

    fake_response = MagicMock()
    fake_response.structures = []
    fake_response.confidence_scores = []
    fake_response.ptm_scores = []

    captured_request = {}

    async def fake_predict(request, **kwargs):
        captured_request["request"] = request
        return fake_response

    with patch("folding.boltz.Boltz2Client") as mock_client_cls:
        mock_client = mock_client_cls.return_value
        mock_client.predict = AsyncMock(side_effect=fake_predict)

        predict_structure(
            "MKV", tmp_path / "out", "test_protein",
            templates=[StructureTemplate(structure_path=template_path, chain_id="A")],
        )

    polymer = captured_request["request"].polymers[0]
    assert len(polymer.structural_templates) == 1
    assert polymer.structural_templates[0].structure == "data_template\n#\n"
    assert polymer.structural_templates[0].chain_id == "A"


def test_predict_structure_rejects_more_than_four_templates(tmp_path, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "nvapi-test")
    template_path = tmp_path / "template.cif"
    template_path.write_text("data_template\n#\n")

    with pytest.raises(ValueError):
        predict_structure(
            "MKV", tmp_path / "out", "test_protein",
            templates=[StructureTemplate(structure_path=template_path)] * 5,
        )

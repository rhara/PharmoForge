import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from proteinanalyzer.cli import protein_info_cmd


def test_protein_info_accepts_accession_directly(tmp_path):
    output = tmp_path / "egfr.json"
    runner = CliRunner()
    with (
        patch("proteinanalyzer.cli.fetch_protein_info", return_value={"accession": "P00533"}) as mock_fetch,
        patch("proteinanalyzer.cli.entry_name_to_accession") as mock_resolve,
        patch("proteinanalyzer.cli.write_protein_info_json") as mock_write,
    ):
        result = runner.invoke(protein_info_cmd, ["P00533", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_not_called()
    mock_fetch.assert_called_once_with("P00533")
    mock_write.assert_called_once_with({"accession": "P00533"}, output)


def test_protein_info_resolves_entry_name_to_accession(tmp_path):
    output = tmp_path / "egfr.json"
    runner = CliRunner()
    with (
        patch("proteinanalyzer.cli.entry_name_to_accession", return_value="P00533") as mock_resolve,
        patch("proteinanalyzer.cli.fetch_protein_info", return_value={"accession": "P00533"}) as mock_fetch,
        patch("proteinanalyzer.cli.write_protein_info_json") as mock_write,
    ):
        result = runner.invoke(protein_info_cmd, ["EGFR_HUMAN", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_called_once_with("EGFR_HUMAN")
    mock_fetch.assert_called_once_with("P00533")
    mock_write.assert_called_once_with({"accession": "P00533"}, output)


def test_protein_info_prints_json_to_stdout_without_output():
    runner = CliRunner()
    with (
        patch("proteinanalyzer.cli.fetch_protein_info", return_value={"accession": "P00533"}) as mock_fetch,
        patch("proteinanalyzer.cli.write_protein_info_json") as mock_write,
    ):
        result = runner.invoke(protein_info_cmd, ["P00533"])

    assert result.exit_code == 0, result.output
    mock_fetch.assert_called_once_with("P00533")
    mock_write.assert_not_called()
    assert json.loads(result.output) == {"accession": "P00533"}

from unittest.mock import patch

from click.testing import CliRunner

from fetcher.cli import fetch_cmd


def test_fetch_activities_dispatch(tmp_path):
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_chembl_target_id", return_value="CHEMBL331") as mock_resolve,
        patch("fetcher.cli.chembl.fetch_activities", return_value=[{"molecule_chembl_id": "CHEMBL1"}]) as mock_fetch,
        patch("fetcher.cli.chembl.write_activities_tsv") as mock_write,
    ):
        result = runner.invoke(fetch_cmd, ["activities=CDK4_HUMAN", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_called_once_with("CDK4_HUMAN")
    mock_fetch.assert_called_once_with("CHEMBL331")
    mock_write.assert_called_once()


def test_fetch_structure_dispatch(tmp_path):
    output = tmp_path / "9csk.cif"
    runner = CliRunner()
    with patch("fetcher.cli.structures.fetch_structure") as mock_fetch_structure:
        result = runner.invoke(fetch_cmd, ["structure=9CSK", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_fetch_structure.assert_called_once_with("9CSK", output)


def test_invalid_spec_missing_equals():
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["activities", "--output", "out.tsv"])
    assert result.exit_code != 0


def test_unsupported_data_type():
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["compounds=CDK4_HUMAN", "--output", "out.tsv"])
    assert result.exit_code != 0

from unittest.mock import patch

from click.testing import CliRunner

from fetcher.cli import fetch_cmd


def test_fetch_activities_dispatch(tmp_path):
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_chembl_target_id", return_value="CHEMBL331") as mock_resolve,
        patch("fetcher.cli.fetch_activities", return_value=[{"molecule_chembl_id": "CHEMBL1"}]) as mock_fetch,
        patch(
            "fetcher.cli.activities.standardize_and_aggregate",
            return_value=[{"smiles": "CCO"}],
        ) as mock_aggregate,
        patch("fetcher.cli.activities.write_activities_tsv") as mock_write,
    ):
        result = runner.invoke(fetch_cmd, ["activities=CDK4_HUMAN", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_called_once_with("CDK4_HUMAN")
    mock_fetch.assert_called_once_with("CHEMBL331")
    mock_aggregate.assert_called_once_with([{"molecule_chembl_id": "CHEMBL1"}])
    mock_write.assert_called_once_with([{"smiles": "CCO"}], output)


def test_fetch_activities_rejects_af_flag(tmp_path):
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["activities=CDK4_HUMAN", "--af", "--output", str(output)])
    assert result.exit_code != 0
    assert "--af" in result.output


def test_fetch_structure_dispatch_infers_format(tmp_path):
    output = tmp_path / "9csk.cif"
    runner = CliRunner()
    with patch("fetcher.cli.fetch_structure") as mock_fetch_structure:
        result = runner.invoke(fetch_cmd, ["structure=9CSK", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_fetch_structure.assert_called_once_with("9CSK", output, fmt=None)


def test_fetch_structure_dispatch_explicit_type(tmp_path):
    output = tmp_path / "9csk.cif"
    runner = CliRunner()
    with patch("fetcher.cli.fetch_structure") as mock_fetch_structure:
        result = runner.invoke(
            fetch_cmd, ["structure=9CSK", "--type=pdb", "--output", str(output)]
        )

    assert result.exit_code == 0, result.output
    mock_fetch_structure.assert_called_once_with("9CSK", output, fmt="pdb")


def test_fetch_structures_dispatch(tmp_path):
    output_dir = tmp_path / "data"
    runner = CliRunner()
    with patch("fetcher.cli.fetch_structures") as mock_fetch_structures:
        result = runner.invoke(
            fetch_cmd,
            ["structures=6P8F,7SJ3,9CSK", "--type", "pdb", "--output", str(output_dir)],
        )

    assert result.exit_code == 0, result.output
    mock_fetch_structures.assert_called_once_with(["6P8F", "7SJ3", "9CSK"], output_dir, "pdb")


def test_fetch_structures_requires_type(tmp_path):
    output_dir = tmp_path / "data"
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["structures=6P8F,7SJ3", "--output", str(output_dir)])

    assert result.exit_code != 0
    assert "--type" in result.output


def test_fetch_structure_af_dispatch_infers_format(tmp_path):
    output = tmp_path / "p61626.cif"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_uniprot_accession", return_value="P61626") as mock_resolve,
        patch("fetcher.cli.fetch_af_structure") as mock_fetch,
    ):
        result = runner.invoke(fetch_cmd, ["structure=P61626", "--af", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_called_once_with("P61626")
    mock_fetch.assert_called_once_with("P61626", output, fmt=None)


def test_fetch_structures_af_dispatch(tmp_path):
    output_dir = tmp_path / "data"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_uniprot_accession", side_effect=["P61626", "P11802"]) as mock_resolve,
        patch("fetcher.cli.fetch_af_structures") as mock_fetch,
    ):
        result = runner.invoke(
            fetch_cmd,
            ["structures=P61626,P11802", "--af", "--type", "pdb", "--output", str(output_dir)],
        )

    assert result.exit_code == 0, result.output
    assert mock_resolve.call_args_list == [(("P61626",),), (("P11802",),)]
    mock_fetch.assert_called_once_with(["P61626", "P11802"], output_dir, "pdb")


def test_fetch_structure_af_dispatch_resolves_entry_name(tmp_path):
    output = tmp_path / "tyk2.cif"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_uniprot_accession", return_value="P29597") as mock_resolve,
        patch("fetcher.cli.fetch_af_structure") as mock_fetch,
    ):
        result = runner.invoke(fetch_cmd, ["structure=TYK2_HUMAN", "--af", "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_called_once_with("TYK2_HUMAN")
    mock_fetch.assert_called_once_with("P29597", output, fmt=None)


def test_fetch_structures_af_dispatch_resolves_mixed_identifiers(tmp_path):
    output_dir = tmp_path / "data"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_uniprot_accession", side_effect=["P61626", "P11802"]) as mock_resolve,
        patch("fetcher.cli.fetch_af_structures") as mock_fetch,
    ):
        result = runner.invoke(
            fetch_cmd,
            ["structures=P61626,CDK4_HUMAN", "--af", "--type", "pdb", "--output", str(output_dir)],
        )

    assert result.exit_code == 0, result.output
    assert mock_resolve.call_args_list == [(("P61626",),), (("CDK4_HUMAN",),)]
    mock_fetch.assert_called_once_with(["P61626", "P11802"], output_dir, "pdb")


def test_fetch_structures_af_requires_type(tmp_path):
    output_dir = tmp_path / "data"
    runner = CliRunner()
    result = runner.invoke(
        fetch_cmd, ["structures=P61626,P11802", "--af", "--output", str(output_dir)]
    )
    assert result.exit_code != 0
    assert "--type" in result.output


def test_invalid_spec_missing_equals():
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["activities", "--output", "out.tsv"])
    assert result.exit_code != 0


def test_unsupported_data_type():
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["compounds=CDK4_HUMAN", "--output", "out.tsv"])
    assert result.exit_code != 0

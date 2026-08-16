from unittest.mock import patch

from click.testing import CliRunner

from fetcher.cli import fetch_cmd


def test_fetch_activity_dispatch(tmp_path):
    outdir = tmp_path / "data"
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
        result = runner.invoke(fetch_cmd, ["activity=CDK4_HUMAN", "--outdir", str(outdir)])

    assert result.exit_code == 0, result.output
    mock_resolve.assert_called_once_with("CDK4_HUMAN")
    mock_fetch.assert_called_once_with("CHEMBL331")
    mock_aggregate.assert_called_once_with([{"molecule_chembl_id": "CHEMBL1"}])
    mock_write.assert_called_once_with([{"smiles": "CCO"}], outdir / "CDK4_HUMAN_activity.tsv")


def test_fetch_activity_rejects_af_flag(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["activity=CDK4_HUMAN", "--af", "--outdir", str(outdir)])
    assert result.exit_code != 0
    assert "--af" in result.output


def test_fetch_structure_single_dispatch_default_format(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    with patch("fetcher.cli.fetch_structure") as mock_fetch_structure:
        result = runner.invoke(fetch_cmd, ["structure=9CSK", "--outdir", str(outdir)])

    assert result.exit_code == 0, result.output
    mock_fetch_structure.assert_called_once_with("9CSK", outdir / "9CSK.cif", fmt="cif")


def test_fetch_structure_single_dispatch_explicit_type(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    with patch("fetcher.cli.fetch_structure") as mock_fetch_structure:
        result = runner.invoke(
            fetch_cmd, ["structure=9csk", "--type=pdb", "--outdir", str(outdir)]
        )

    assert result.exit_code == 0, result.output
    mock_fetch_structure.assert_called_once_with("9csk", outdir / "9CSK.pdb", fmt="pdb")


def test_fetch_structure_multiple_dispatch(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    with patch("fetcher.cli.fetch_structure") as mock_fetch_structure:
        result = runner.invoke(
            fetch_cmd,
            ["structure=6P8F,7SJ3,9CSK", "--type", "pdb", "--outdir", str(outdir)],
        )

    assert result.exit_code == 0, result.output
    assert mock_fetch_structure.call_args_list == [
        (("6P8F", outdir / "6P8F.pdb"), {"fmt": "pdb"}),
        (("7SJ3", outdir / "7SJ3.pdb"), {"fmt": "pdb"}),
        (("9CSK", outdir / "9CSK.pdb"), {"fmt": "pdb"}),
    ]


def test_fetch_structure_af_single_dispatch(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_uniprot_accession", return_value="P61626") as mock_resolve_acc,
        patch("fetcher.cli.fetch_af_structure") as mock_fetch,
    ):
        result = runner.invoke(fetch_cmd, ["structure=P61626", "--af", "--outdir", str(outdir)])

    assert result.exit_code == 0, result.output
    mock_resolve_acc.assert_called_once_with("P61626")
    mock_fetch.assert_called_once_with("P61626", outdir / "P61626_AF.cif", fmt="cif")


def test_fetch_structure_af_dispatch_resolves_accession_from_entry_name(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    with (
        patch("fetcher.cli.resolve_uniprot_accession", return_value="P29597") as mock_resolve_acc,
        patch("fetcher.cli.fetch_af_structure") as mock_fetch,
    ):
        result = runner.invoke(
            fetch_cmd, ["structure=TYK2_HUMAN", "--af", "--type=cif", "--outdir", str(outdir)]
        )

    assert result.exit_code == 0, result.output
    mock_resolve_acc.assert_called_once_with("TYK2_HUMAN")
    mock_fetch.assert_called_once_with("P29597", outdir / "P29597_AF.cif", fmt="cif")


def test_fetch_structure_af_multiple_dispatch_resolves_mixed_identifiers(tmp_path):
    outdir = tmp_path / "data"
    runner = CliRunner()
    with (
        patch(
            "fetcher.cli.resolve_uniprot_accession", side_effect=["P61626", "P11802"]
        ) as mock_resolve_acc,
        patch("fetcher.cli.fetch_af_structure") as mock_fetch,
    ):
        result = runner.invoke(
            fetch_cmd,
            ["structure=P61626,CDK4_HUMAN", "--af", "--type", "pdb", "--outdir", str(outdir)],
        )

    assert result.exit_code == 0, result.output
    assert mock_resolve_acc.call_args_list == [(("P61626",),), (("CDK4_HUMAN",),)]
    assert mock_fetch.call_args_list == [
        (("P61626", outdir / "P61626_AF.pdb"), {"fmt": "pdb"}),
        (("P11802", outdir / "P11802_AF.pdb"), {"fmt": "pdb"}),
    ]


def test_outdir_is_required():
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["structure=9CSK"])
    assert result.exit_code != 0
    assert "--outdir" in result.output or "-o" in result.output


def test_invalid_spec_missing_equals(tmp_path):
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["activity", "--outdir", str(tmp_path)])
    assert result.exit_code != 0


def test_unsupported_data_type(tmp_path):
    runner = CliRunner()
    result = runner.invoke(fetch_cmd, ["compounds=CDK4_HUMAN", "--outdir", str(tmp_path)])
    assert result.exit_code != 0

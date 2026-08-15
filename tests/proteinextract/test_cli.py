from unittest.mock import patch

from click.testing import CliRunner

from proteinextract.cli import protein_extract_cmd


def test_protein_extract_parses_chains_and_passes_options(tmp_path):
    input_path = tmp_path / "input.pdb"
    input_path.write_text("ATOM\n")
    output = tmp_path / "out.cif"
    runner = CliRunner()
    with patch("proteinextract.cli.extract_structure") as mock_extract:
        result = runner.invoke(
            protein_extract_cmd,
            [str(input_path), "--chains=A,G,H,I,J", "--remove-water", "--output", str(output)],
        )

    assert result.exit_code == 0, result.output
    mock_extract.assert_called_once_with(
        input_path, output, chains=["A", "G", "H", "I", "J"], remove_water=True
    )


def test_protein_extract_defaults_no_chains_no_remove_water(tmp_path):
    input_path = tmp_path / "input.pdb"
    input_path.write_text("ATOM\n")
    output = tmp_path / "out.pdb"
    runner = CliRunner()
    with patch("proteinextract.cli.extract_structure") as mock_extract:
        result = runner.invoke(protein_extract_cmd, [str(input_path), "--output", str(output)])

    assert result.exit_code == 0, result.output
    mock_extract.assert_called_once_with(input_path, output, chains=None, remove_water=False)


def test_protein_extract_requires_existing_input():
    runner = CliRunner()
    result = runner.invoke(protein_extract_cmd, ["nonexistent.pdb", "--output", "out.pdb"])
    assert result.exit_code != 0

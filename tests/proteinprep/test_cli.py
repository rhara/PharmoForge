from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from proteinprep.cli import prep_protein_cmd


def test_prep_protein_dock_mode_passes_ph_none(tmp_path):
    input_path = tmp_path / "input.pdb"
    input_path.write_text("ATOM\n")
    output = tmp_path / "out.pdb"
    runner = CliRunner()
    with patch("proteinprep.cli.repair_structure") as mock_repair:
        result = runner.invoke(
            prep_protein_cmd, [str(input_path), "--output", str(output)]
        )

    assert result.exit_code == 0, result.output
    mock_repair.assert_called_once_with(input_path, output, ph=None)


def test_prep_protein_md_mode_passes_ph(tmp_path):
    input_path = tmp_path / "input.pdb"
    input_path.write_text("ATOM\n")
    output = tmp_path / "out.pdb"
    runner = CliRunner()
    with patch("proteinprep.cli.repair_structure") as mock_repair:
        result = runner.invoke(
            prep_protein_cmd,
            [str(input_path), "--output", str(output), "--mode", "md", "--ph", "7.4"],
        )

    assert result.exit_code == 0, result.output
    mock_repair.assert_called_once_with(input_path, output, ph=7.4)


def test_prep_protein_requires_existing_input():
    runner = CliRunner()
    result = runner.invoke(
        prep_protein_cmd, ["nonexistent.pdb", "--output", "out.pdb"]
    )
    assert result.exit_code != 0

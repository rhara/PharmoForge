from unittest.mock import patch

from click.testing import CliRunner

from alignview.cli import align_view_cmd


def test_align_view_passes_paths_and_options(tmp_path):
    p1 = tmp_path / "a.pdb"
    p2 = tmp_path / "b.cif"
    p1.write_text("ATOM\n")
    p2.write_text("data_b\n")
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(
            align_view_cmd,
            [
                str(p1),
                str(p2),
                "--method",
                "cealign",
                "--align-margin",
                "30",
                "--pymol-env",
                "custom-pymol",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with(
        [p1, p2], method="cealign", pymol_env="custom-pymol", align_margin=30
    )


def test_align_view_defaults(tmp_path):
    p1 = tmp_path / "a.pdb"
    p1.write_text("ATOM\n")
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(align_view_cmd, [str(p1)])

    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with([p1], method="align", pymol_env="pymol", align_margin=20)


def test_align_view_accepts_number_method(tmp_path):
    p1 = tmp_path / "a.pdb"
    p2 = tmp_path / "b.pdb"
    p1.write_text("ATOM\n")
    p2.write_text("ATOM\n")
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(align_view_cmd, [str(p1), str(p2), "--method", "number"])

    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with([p1, p2], method="number", pymol_env="pymol", align_margin=20)


def test_align_view_rejects_unknown_method(tmp_path):
    p1 = tmp_path / "a.pdb"
    p1.write_text("ATOM\n")
    runner = CliRunner()
    result = runner.invoke(align_view_cmd, [str(p1), "--method", "bogus"])
    assert result.exit_code != 0


def test_align_view_requires_existing_paths():
    runner = CliRunner()
    result = runner.invoke(align_view_cmd, ["nonexistent.pdb"])
    assert result.exit_code != 0


def test_align_view_requires_at_least_one_path():
    runner = CliRunner()
    result = runner.invoke(align_view_cmd, [])
    assert result.exit_code != 0

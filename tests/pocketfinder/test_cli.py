from unittest.mock import patch

from click.testing import CliRunner

from pocket import Pocket
from pocketfinder.cli import find_pocket_cmd, show_pocket_cmd, view_pocket_cmd


def _make_pockets(n: int) -> list[Pocket]:
    return [
        Pocket(
            pocket_id=i,
            score=1.0 / i,
            druggability_score=0.5,
            n_alpha_spheres=10,
            volume=100.0,
            total_sasa=10.0,
            polar_sasa=5.0,
            apolar_sasa=5.0,
            hydrophobicity_score=1.0,
            residues=[],
        )
        for i in range(1, n + 1)
    ]


def runner_invoke(structure_path, output_dir, extra=None):
    runner = CliRunner()
    args = [str(structure_path), "--outdir", str(output_dir)] + (extra or [])
    return runner.invoke(find_pocket_cmd, args)


def test_find_pocket_writes_all_pockets_by_default(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")
    output_dir = tmp_path / "out"
    pockets = _make_pockets(3)

    with (
        patch("pocketfinder.cli.run_fpocket", return_value=pockets) as mock_run,
        patch("pocketfinder.cli.write_pockets_json") as mock_write,
    ):
        result = runner_invoke(structure_path, output_dir)

    assert result.exit_code == 0, result.output
    mock_run.assert_called_once_with(structure_path, output_dir)
    mock_write.assert_called_once_with("input.pdb", pockets, output_dir / "pockets.json")


def test_find_pocket_applies_top_n(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")
    output_dir = tmp_path / "out"
    pockets = _make_pockets(5)

    with (
        patch("pocketfinder.cli.run_fpocket", return_value=pockets),
        patch("pocketfinder.cli.write_pockets_json") as mock_write,
    ):
        result = runner_invoke(structure_path, output_dir, extra=["--top", "2"])

    assert result.exit_code == 0, result.output
    mock_write.assert_called_once_with("input.pdb", pockets[:2], output_dir / "pockets.json")


def test_find_pocket_requires_existing_structure():
    runner = CliRunner()
    result = runner.invoke(find_pocket_cmd, ["nonexistent.pdb", "--outdir", "out"])
    assert result.exit_code != 0


def test_view_pocket_reads_json_and_launches_view(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")
    pockets_path = tmp_path / "pockets.json"
    pockets_path.write_text("{}")
    pockets = _make_pockets(3)

    runner = CliRunner()
    with (
        patch("pocketfinder.cli.read_pockets_json", return_value=pockets) as mock_read,
        patch("pocketfinder.cli.launch_pocket_view") as mock_launch,
    ):
        result = runner.invoke(
            view_pocket_cmd,
            [str(structure_path), "--pockets", str(pockets_path), "--top", "2", "--pymol-env", "custom"],
        )

    assert result.exit_code == 0, result.output
    mock_read.assert_called_once_with(pockets_path)
    mock_launch.assert_called_once_with(structure_path, pockets, pymol_env="custom", top_n=2)


def test_view_pocket_requires_existing_pockets_file(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")

    runner = CliRunner()
    result = runner.invoke(
        view_pocket_cmd, [str(structure_path), "--pockets", str(tmp_path / "nonexistent.json")]
    )
    assert result.exit_code != 0


def test_show_pocket_prints_table_to_stdout_by_default(tmp_path):
    pockets_path = tmp_path / "pockets.json"
    pockets_path.write_text("{}")
    pockets = _make_pockets(2)

    runner = CliRunner()
    with (
        patch("pocketfinder.cli.read_pockets_json", return_value=pockets) as mock_read,
        patch("pocketfinder.cli.write_pockets_table") as mock_write,
    ):
        result = runner.invoke(show_pocket_cmd, [str(pockets_path)])

    assert result.exit_code == 0, result.output
    mock_read.assert_called_once_with(pockets_path)
    mock_write.assert_not_called()
    assert "pocket_id" in result.output
    assert "1\t" in result.output


def test_show_pocket_writes_table_when_output_given(tmp_path):
    pockets_path = tmp_path / "pockets.json"
    pockets_path.write_text("{}")
    output_path = tmp_path / "pockets.tsv"
    pockets = _make_pockets(2)

    runner = CliRunner()
    with (
        patch("pocketfinder.cli.read_pockets_json", return_value=pockets),
        patch("pocketfinder.cli.write_pockets_table") as mock_write,
    ):
        result = runner.invoke(show_pocket_cmd, [str(pockets_path), "--output", str(output_path)])

    assert result.exit_code == 0, result.output
    mock_write.assert_called_once_with(pockets, output_path)


def test_show_pocket_requires_existing_pockets_file():
    runner = CliRunner()
    result = runner.invoke(show_pocket_cmd, ["nonexistent.json"])
    assert result.exit_code != 0

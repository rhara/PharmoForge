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


def test_align_view_indir_resolves_stems_with_auto_extension(tmp_path):
    (tmp_path / "a.cif").write_text("data_a\n")
    (tmp_path / "b.pdb").write_text("ATOM\n")
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(align_view_cmd, ["--indir", str(tmp_path), "a", "b"])

    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with(
        [tmp_path / "a.cif", tmp_path / "b.pdb"], method="align", pymol_env="pymol", align_margin=20
    )


def test_align_view_indir_prefers_cif_over_pdb(tmp_path):
    (tmp_path / "a.cif").write_text("data_a\n")
    (tmp_path / "a.pdb").write_text("ATOM\n")
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(align_view_cmd, ["--indir", str(tmp_path), "a"])

    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with([tmp_path / "a.cif"], method="align", pymol_env="pymol", align_margin=20)


def test_align_view_indir_can_be_repeated(tmp_path):
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "a.cif").write_text("data_a\n")
    (dir2 / "b.pdb").write_text("ATOM\n")
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(
            align_view_cmd, ["--indir", str(dir1), "a", "--indir", str(dir2), "b"]
        )

    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with(
        [dir1 / "a.cif", dir2 / "b.pdb"], method="align", pymol_env="pymol", align_margin=20
    )


def test_align_view_slash_containing_token_ignores_indir(tmp_path):
    outside_dir = tmp_path / "outside"
    outside_dir.mkdir()
    outside_file = outside_dir / "c.cif"
    outside_file.write_text("data_c\n")
    indir = tmp_path / "indir"
    indir.mkdir()
    runner = CliRunner()
    with patch("alignview.cli.launch_alignment_view") as mock_launch:
        result = runner.invoke(
            align_view_cmd, ["--indir", str(indir), str(outside_file)]
        )

    # 絶対パス("/"を含む)指定はindirを無視してそのまま使われる
    assert result.exit_code == 0, result.output
    mock_launch.assert_called_once_with(
        [outside_file], method="align", pymol_env="pymol", align_margin=20
    )


def test_align_view_indir_missing_directory_value():
    runner = CliRunner()
    result = runner.invoke(align_view_cmd, ["--indir"])
    assert result.exit_code != 0


def test_align_view_indir_stem_not_found(tmp_path):
    runner = CliRunner()
    result = runner.invoke(align_view_cmd, ["--indir", str(tmp_path), "missing"])
    assert result.exit_code != 0

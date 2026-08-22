from unittest.mock import patch

import pytest

from docking.vina import calc_search_box, parse_vina_output, run_vina

_SAMPLE_OUTPUT_PDBQT = """\
MODEL 1
REMARK VINA RESULT:    -7.3      0.000      0.000
ATOM      1  C   LIG A   1       1.000   1.000   1.000  1.00  0.00     0.000 C
ENDMDL
MODEL 2
REMARK VINA RESULT:    -6.8      1.421      2.103
ATOM      1  C   LIG A   1       1.100   1.200   1.300  1.00  0.00     0.000 C
ENDMDL
"""


def test_calc_search_box_returns_center_and_padded_size():
    coords = [(0.0, 0.0, 0.0), (10.0, 4.0, 2.0)]

    center, size = calc_search_box(coords, padding=5.0)

    assert center == (5.0, 2.0, 1.0)
    assert size == (10.0 + 10.0, 4.0 + 10.0, 2.0 + 10.0)


def test_run_vina_invokes_mamba_run_vina_and_parses_poses(tmp_path):
    output_path = tmp_path / "docked.pdbqt"

    def _fake_run(command, capture_output, text):
        output_path.write_text(_SAMPLE_OUTPUT_PDBQT)

        class _Result:
            returncode = 0
            stderr = ""

        return _Result()

    with (
        patch("docking.vina.subprocess.run", side_effect=_fake_run) as mock_run,
        patch("docking.vina.shutil.which", return_value="/usr/bin/mamba"),
    ):
        result = run_vina(
            rigid_pdbqt=tmp_path / "rigid.pdbqt",
            ligand_pdbqt=tmp_path / "ligand.pdbqt",
            center=(1.0, 2.0, 3.0),
            size=(20.0, 20.0, 20.0),
            output_path=output_path,
            flex_pdbqt=tmp_path / "flex.pdbqt",
            vina_env="vina",
        )

    command = mock_run.call_args.args[0]
    assert command[:4] == ["/usr/bin/mamba", "run", "-n", "vina"]
    assert command[4] == "vina"
    assert "--flex" in command

    assert [p.affinity for p in result.poses] == [-7.3, -6.8]
    assert result.best_affinity == -7.3
    assert result.poses[1].rmsd_lb == 1.421


def test_run_vina_raises_when_mamba_not_found(tmp_path):
    with patch("docking.vina.shutil.which", return_value=None):
        with pytest.raises(RuntimeError):
            run_vina(
                rigid_pdbqt=tmp_path / "rigid.pdbqt",
                ligand_pdbqt=tmp_path / "ligand.pdbqt",
                center=(0.0, 0.0, 0.0),
                size=(20.0, 20.0, 20.0),
                output_path=tmp_path / "docked.pdbqt",
            )


def test_parse_vina_output_reads_cached_file_without_rerunning(tmp_path):
    output_path = tmp_path / "docked.pdbqt"
    output_path.write_text(_SAMPLE_OUTPUT_PDBQT)

    poses = parse_vina_output(output_path)

    assert [p.affinity for p in poses] == [-7.3, -6.8]
    assert poses[0].mode == 1


def test_parse_vina_output_raises_when_no_poses_found(tmp_path):
    output_path = tmp_path / "empty.pdbqt"
    output_path.write_text("MODEL 1\nENDMDL\n")

    with pytest.raises(ValueError):
        parse_vina_output(output_path)


def test_run_vina_raises_on_nonzero_exit(tmp_path):
    def _fake_run(command, capture_output, text):
        class _Result:
            returncode = 1
            stderr = "boom"

        return _Result()

    with (
        patch("docking.vina.subprocess.run", side_effect=_fake_run),
        patch("docking.vina.shutil.which", return_value="/usr/bin/mamba"),
    ):
        with pytest.raises(RuntimeError):
            run_vina(
                rigid_pdbqt=tmp_path / "rigid.pdbqt",
                ligand_pdbqt=tmp_path / "ligand.pdbqt",
                center=(0.0, 0.0, 0.0),
                size=(20.0, 20.0, 20.0),
                output_path=tmp_path / "docked.pdbqt",
            )

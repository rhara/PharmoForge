from unittest.mock import patch

import numpy as np

from alignview.view import build_pymol_script, launch_alignment_view
from structfit import FitResult


def test_build_pymol_script_align_embeds_range_restricted_python_block(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p2 = tmp_path / "mobile.pdb"
    p1.write_text("")
    p2.write_text("")

    script = build_pymol_script([p1, p2], method="align", align_margin=20)

    assert f"load {p1.resolve()}, ref" in script
    assert f"load {p2.resolve()}, mobile" in script
    assert 'cmd.iterate("mobile and polymer.protein and name CA"' in script
    assert "_lo, _hi = min(_resis) - 20, max(_resis) + 20" in script
    assert 'cmd.align("mobile", _target_sel)' in script
    assert 'cmd.align("mobile", "ref")' in script  # フォールバック時の呼び出し
    assert "zoom all" in script


def test_build_pymol_script_super_and_cealign_are_plain_commands(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p2 = tmp_path / "mobile.pdb"
    p1.write_text("")
    p2.write_text("")

    super_script = build_pymol_script([p1, p2], method="super")
    assert "super mobile, ref" in super_script
    assert "python" not in super_script

    cealign_script = build_pymol_script([p1, p2], method="cealign")
    assert "cealign ref, mobile" in cealign_script
    assert "python" not in cealign_script


def test_build_pymol_script_dedups_duplicate_stems(tmp_path):
    d1 = tmp_path / "a"
    d2 = tmp_path / "b"
    d1.mkdir()
    d2.mkdir()
    p1 = d1 / "same.pdb"
    p2 = d2 / "same.pdb"
    p1.write_text("")
    p2.write_text("")

    script = build_pymol_script([p1, p2], method="cealign")

    assert f"load {p1.resolve()}, same" in script
    assert "same_2" in script
    assert "cealign same, same_2" in script


def test_build_pymol_script_number_embeds_precomputed_transform(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p2 = tmp_path / "mobile.pdb"
    p1.write_text("")
    p2.write_text("")

    fake_result = FitResult(
        matrix=np.eye(4),
        rmsd=0.5,
        n_residues=100,
        mobile_chain="A",
        target_chain="B",
    )
    with patch("alignview.view.fit_by_residue_number", return_value=fake_result) as mock_fit:
        script = build_pymol_script([p1, p2], method="number")

    mock_fit.assert_called_once_with(p2, p1)
    assert 'cmd.transform_object("mobile",' in script
    assert "transpose=0" in script
    assert "RMSD=%.3f" in script
    assert "zoom all" in script


def test_build_pymol_script_number_falls_back_on_fit_failure(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p2 = tmp_path / "mobile.pdb"
    p1.write_text("")
    p2.write_text("")

    with patch("alignview.view.fit_by_residue_number", side_effect=ValueError("no common residues")):
        script = build_pymol_script([p1, p2], method="number")

    assert "transform_object" not in script
    assert "number-based fit failed for mobile" in script
    assert "zoom all" in script


def test_build_pymol_script_single_structure_has_no_align(tmp_path):
    p1 = tmp_path / "only.pdb"
    p1.write_text("")

    script = build_pymol_script([p1], method="align")

    assert "load" in script
    assert "align" not in script
    assert "zoom all" in script


def test_launch_alignment_view_builds_script_and_runs_pymol(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p1.write_text("")

    with patch("alignview.view.run_pymol_script") as mock_run:
        launch_alignment_view([p1], method="align", pymol_env="custom-pymol")

    args, kwargs = mock_run.call_args
    assert f"load {p1.resolve()}, ref" in args[0]
    assert kwargs == {"pymol_env": "custom-pymol"}

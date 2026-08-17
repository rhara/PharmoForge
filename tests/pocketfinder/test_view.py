from unittest.mock import patch

from pocket import Pocket, PocketResidue
from pocketfinder.view import build_pocket_view_script, launch_pocket_view


def _pocket(pocket_id: int, residues: list[PocketResidue]) -> Pocket:
    return Pocket(
        pocket_id=pocket_id,
        score=1.0 / pocket_id,
        druggability_score=0.5,
        n_alpha_spheres=10,
        volume=100.0,
        total_sasa=10.0,
        polar_sasa=5.0,
        apolar_sasa=5.0,
        hydrophobicity_score=1.0,
        residues=residues,
    )


def test_build_pocket_view_script_selects_residues_by_chain(tmp_path):
    structure_path = tmp_path / "P24941_AF.cif"
    structure_path.write_text("data_x\n")
    pockets = [
        _pocket(1, [PocketResidue("A", 14, "THR"), PocketResidue("A", 15, "TYR")]),
        _pocket(2, [PocketResidue("B", 3, "VAL")]),
    ]

    script = build_pocket_view_script(structure_path, pockets)

    assert f"load {structure_path.resolve()}, P24941_AF" in script
    assert "select pocket_1, P24941_AF and ((chain A and resi 14+15))" in script
    assert "select pocket_2, P24941_AF and ((chain B and resi 3))" in script
    assert "zoom (pocket_1 or pocket_2)" in script


def test_build_pocket_view_script_applies_top_n(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")
    pockets = [
        _pocket(1, [PocketResidue("A", 1, "GLY")]),
        _pocket(2, [PocketResidue("A", 2, "GLY")]),
        _pocket(3, [PocketResidue("A", 3, "GLY")]),
    ]

    script = build_pocket_view_script(structure_path, pockets, top_n=1)

    assert "pocket_1" in script
    assert "pocket_2" not in script
    assert "pocket_3" not in script


def test_build_pocket_view_script_skips_pockets_without_residues(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")
    pockets = [_pocket(1, [])]

    script = build_pocket_view_script(structure_path, pockets)

    assert "pocket_1" not in script
    assert "zoom input" in script


def test_launch_pocket_view_builds_script_and_runs_pymol(tmp_path):
    structure_path = tmp_path / "input.pdb"
    structure_path.write_text("ATOM\n")
    pockets = [_pocket(1, [PocketResidue("A", 1, "GLY")])]

    with patch("pocketfinder.view.run_pymol_script") as mock_run:
        launch_pocket_view(structure_path, pockets, pymol_env="custom-pymol", top_n=1)

    args, kwargs = mock_run.call_args
    assert "pocket_1" in args[0]
    assert kwargs == {"pymol_env": "custom-pymol"}

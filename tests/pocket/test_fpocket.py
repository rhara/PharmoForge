from unittest.mock import patch

import pytest

from pocket.fpocket import Pocket, PocketResidue, _extract_residues, _parse_info_file, run_fpocket

_INFO_TEXT = """Pocket 1 :
\tScore : \t0.549
\tDruggability Score : \t0.853
\tNumber of Alpha Spheres : \t55
\tTotal SASA : \t64.374
\tPolar SASA : \t28.146
\tApolar SASA : \t36.229
\tVolume : \t441.369
\tMean local hydrophobic density : \t33.737
\tMean alpha sphere radius :\t4.049
\tHydrophobicity score:\t46.429
\tVolume score: \t 4.500
\tCharge score :\t -2

Pocket 2 :
\tScore : \t0.255
\tDruggability Score : \t0.001
\tNumber of Alpha Spheres : \t17
\tTotal SASA : \t47.177
\tPolar SASA : \t29.062
\tApolar SASA : \t18.114
\tVolume : \t308.239
\tHydrophobicity score:\t41.429
"""


def _atom_line(serial: int, resseq: int, chain: str, resname: str, atom: str = "CA") -> str:
    return (
        f"ATOM  {serial:>5} {atom:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {atom[0]:>2}"
    )


def test_parse_info_file_extracts_known_fields_only(tmp_path):
    info_path = tmp_path / "test_info.txt"
    info_path.write_text(_INFO_TEXT)

    infos = _parse_info_file(info_path)

    assert set(infos) == {1, 2}
    assert infos[1] == {
        "score": 0.549,
        "druggability_score": 0.853,
        "n_alpha_spheres": 55.0,
        "total_sasa": 64.374,
        "polar_sasa": 28.146,
        "apolar_sasa": 36.229,
        "volume": 441.369,
        "hydrophobicity_score": 46.429,
    }
    assert infos[2]["score"] == 0.255


def test_extract_residues_deduplicates_and_sorts(tmp_path):
    atm_path = tmp_path / "pocket1_atm.pdb"
    lines = [
        _atom_line(1, 20, "A", "LEU", atom="CA"),
        _atom_line(2, 20, "A", "LEU", atom="CB"),
        _atom_line(3, 5, "A", "GLY", atom="CA"),
        _atom_line(4, 3, "B", "VAL", atom="CA"),
        "TER",
        "END",
    ]
    atm_path.write_text("\n".join(lines) + "\n")

    residues = _extract_residues(atm_path)

    assert residues == [
        PocketResidue(chain_id="A", resnum=5, resname="GLY"),
        PocketResidue(chain_id="A", resnum=20, resname="LEU"),
        PocketResidue(chain_id="B", resnum=3, resname="VAL"),
    ]


def _write_fpocket_output(work_dir, stem: str):
    out_dir = work_dir / f"{stem}_out"
    (out_dir / "pockets").mkdir(parents=True)
    (out_dir / f"{stem}_info.txt").write_text(_INFO_TEXT)
    (out_dir / "pockets" / "pocket1_atm.pdb").write_text(_atom_line(1, 5, "A", "GLY") + "\n")
    (out_dir / "pockets" / "pocket2_atm.pdb").write_text(_atom_line(1, 9, "A", "SER") + "\n")


def test_run_fpocket_invokes_binary_and_parses_output(tmp_path):
    structure_path = tmp_path / "input" / "test.pdb"
    structure_path.parent.mkdir()
    structure_path.write_text("ATOM\n")
    work_dir = tmp_path / "work"

    def fake_run(*args, **kwargs):
        _write_fpocket_output(work_dir, "test")

        class Result:
            returncode = 0
            stderr = ""

        return Result()

    with patch("pocket.fpocket.subprocess.run", side_effect=fake_run) as mock_run:
        pockets = run_fpocket(structure_path, work_dir)

    args, kwargs = mock_run.call_args
    assert args[0] == ["fpocket", "-f", "test.pdb", "-w", "pdb"]
    assert kwargs["cwd"] == work_dir
    assert (work_dir / "test.pdb").exists()

    assert pockets == [
        Pocket(
            pocket_id=1,
            score=0.549,
            druggability_score=0.853,
            n_alpha_spheres=55,
            volume=441.369,
            total_sasa=64.374,
            polar_sasa=28.146,
            apolar_sasa=36.229,
            hydrophobicity_score=46.429,
            residues=[PocketResidue(chain_id="A", resnum=5, resname="GLY")],
        ),
        Pocket(
            pocket_id=2,
            score=0.255,
            druggability_score=0.001,
            n_alpha_spheres=17,
            volume=308.239,
            total_sasa=47.177,
            polar_sasa=29.062,
            apolar_sasa=18.114,
            hydrophobicity_score=41.429,
            residues=[PocketResidue(chain_id="A", resnum=9, resname="SER")],
        ),
    ]


def test_run_fpocket_raises_on_nonzero_exit(tmp_path):
    structure_path = tmp_path / "test.pdb"
    structure_path.write_text("ATOM\n")

    class Result:
        returncode = 1
        stderr = "Invalid pdb name given."

    with patch("pocket.fpocket.subprocess.run", return_value=Result()):
        with pytest.raises(RuntimeError, match="Invalid pdb name given"):
            run_fpocket(structure_path, tmp_path / "work")

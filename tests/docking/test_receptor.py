import pytest

from docking import prepare_flexible_receptor

# CDK20_HUMAN(AlphaFold予測構造、Q8IZL9)のGLY11-GLU12-GLY13の実座標を抜粋した最小フラグメント。
# GLUは回転可能な側鎖(CB-CG-CD-OE1/OE2)を持つため、フレキシブル化のテストに使う。
FRAGMENT_PDB = """\
ATOM      1  N   GLY A  11      -0.873   3.233  18.545  1.00 76.12           N
ATOM      2  CA  GLY A  11       0.152   3.760  17.654  1.00 76.12           C
ATOM      3  C   GLY A  11       0.506   5.216  17.913  1.00 76.12           C
ATOM      4  O   GLY A  11      -0.144   5.921  18.690  1.00 76.12           O
ATOM      5  N   GLU A  12       1.559   5.653  17.238  1.00 72.38           N
ATOM      6  CA  GLU A  12       1.976   7.045  17.110  1.00 72.38           C
ATOM      7  C   GLU A  12       2.250   7.299  15.631  1.00 72.38           C
ATOM      8  CB  GLU A  12       3.227   7.291  17.961  1.00 72.38           C
ATOM      9  O   GLU A  12       3.081   6.623  15.029  1.00 72.38           O
ATOM     10  CG  GLU A  12       3.456   8.790  18.209  1.00 72.38           C
ATOM     11  CD  GLU A  12       4.594   9.086  19.203  1.00 72.38           C
ATOM     12  OE1 GLU A  12       4.930  10.275  19.343  1.00 72.38           O
ATOM     13  OE2 GLU A  12       5.029   8.148  19.915  1.00 72.38           O
ATOM     14  N   GLY A  13       1.496   8.213  15.030  1.00 63.19           N
ATOM     15  CA  GLY A  13       1.735   8.681  13.671  1.00 63.19           C
ATOM     16  C   GLY A  13       2.440  10.033  13.678  1.00 63.19           C
ATOM     17  O   GLY A  13       2.541  10.693  14.710  1.00 63.19           O
END
"""


def _write_fragment(tmp_path):
    input_path = tmp_path / "fragment.pdb"
    input_path.write_text(FRAGMENT_PDB)
    return input_path


def test_prepare_flexible_receptor_splits_rigid_and_flex(tmp_path):
    input_path = _write_fragment(tmp_path)

    result = prepare_flexible_receptor(input_path, [("A", 12)], tmp_path / "out")

    assert result.rigid_pdbqt.exists()
    assert result.flex_pdbqt is not None
    assert result.flex_pdbqt.exists()
    assert result.polymer_json.exists()
    assert result.n_flexible_residues == 1

    flex_text = result.flex_pdbqt.read_text()
    assert "GLU A  12" in flex_text
    rigid_text = result.rigid_pdbqt.read_text()
    assert "GLY A  11" in rigid_text  # フレキシブル指定していない残基は rigid 側に残る


def test_prepare_flexible_receptor_no_flexible_residues_omits_flex_file(tmp_path):
    input_path = _write_fragment(tmp_path)

    result = prepare_flexible_receptor(input_path, [], tmp_path / "out")

    assert result.flex_pdbqt is None
    assert result.n_flexible_residues == 0


def test_prepare_flexible_receptor_rejects_unknown_residue(tmp_path):
    input_path = _write_fragment(tmp_path)

    with pytest.raises(ValueError):
        prepare_flexible_receptor(input_path, [("A", 999)], tmp_path / "out")

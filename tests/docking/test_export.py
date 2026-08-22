from pathlib import Path

from ligandprep import prepare_ligand_pdbqt

from docking import export_docked_poses, prepare_flexible_receptor

# CDK20_HUMAN(AlphaFold予測構造、Q8IZL9)のGLY11-GLU12-GLY13の実座標を抜粋した最小フラグメント
# (tests/docking/test_receptor.pyと同じ)。GLUは回転可能な側鎖を持つためフレキシブル化のテストに使う。
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


def _build_fake_vina_output(flex_pdbqt_text: str, ligand_pdbqt_text: str, tmp_path) -> Path:
    # 実際のVina出力は1モデル毎にリガンド原子(REMARK VINA RESULT付き)+可動側鎖(BEGIN_RES/END_RES)
    # をMODEL/ENDMDLで囲んで書き出す(docking.vinaのテスト・実データで確認済みの実フォーマット)。
    # ここではreceptorの初期配座をそのまま「ドッキング後」の座標として流用し、1モードだけ組み立てる。
    ligand_body = "\n".join(
        line for line in ligand_pdbqt_text.splitlines() if not line.startswith("REMARK INDEX MAP")
    )
    flex_body = flex_pdbqt_text
    output = (
        "MODEL 1\n"
        "REMARK VINA RESULT:    -5.000      0.000      0.000\n"
        f"{ligand_body}\n"
        f"{flex_body}"
        "ENDMDL\n"
    )
    out_path = tmp_path / "fake_docked.pdbqt"
    out_path.write_text(output)
    return out_path


def test_export_docked_poses_writes_receptor_pdb_and_ligand_sdf(tmp_path):
    input_path = tmp_path / "fragment.pdb"
    input_path.write_text(FRAGMENT_PDB)

    flex_receptor = prepare_flexible_receptor(input_path, [("A", 12)], tmp_path / "receptor")
    ligand_path = prepare_ligand_pdbqt("CO", "methanol", tmp_path / "methanol.pdbqt")

    fake_output = _build_fake_vina_output(
        flex_receptor.flex_pdbqt.read_text(), ligand_path.read_text(), tmp_path
    )

    exported = export_docked_poses(
        polymer_json=flex_receptor.polymer_json,
        vina_output_pdbqt=fake_output,
        output_dir=tmp_path / "exported",
        name="methanol",
    )

    assert len(exported) == 1
    pose = exported[0]
    assert pose.mode == 1
    assert pose.receptor_pdb.exists()
    assert pose.ligand_sdf.exists()

    receptor_text = pose.receptor_pdb.read_text()
    assert "GLY A  11" in receptor_text
    assert "GLU A  12" in receptor_text
    assert "GLY A  13" in receptor_text

    ligand_text = pose.ligand_sdf.read_text()
    assert "V2000" in ligand_text
    assert ligand_text.strip().endswith("$$$$")

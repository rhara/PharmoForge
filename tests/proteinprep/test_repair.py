from proteinprep import repair_structure

# 欠損原子・欠損残基を含む最小限のジペプチド(ALA-GLY)。側鎖・水素は未指定。
MINIMAL_PDB = """\
ATOM      1  N   ALA A   1      -0.966   0.493   1.500  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.257   1.311   1.500  1.00  0.00           C
ATOM      3  C   ALA A   1       1.494   0.431   1.500  1.00  0.00           C
ATOM      4  O   ALA A   1       1.494  -0.802   1.500  1.00  0.00           O
ATOM      5  N   GLY A   2       2.628   1.113   1.500  1.00  0.00           N
ATOM      6  CA  GLY A   2       3.909   0.421   1.500  1.00  0.00           C
ATOM      7  C   GLY A   2       5.100   1.360   1.500  1.00  0.00           C
ATOM      8  O   GLY A   2       5.100   2.593   1.500  1.00  0.00           O
TER       9      GLY A   2
END
"""


def _write_input(tmp_path):
    input_path = tmp_path / "input.pdb"
    input_path.write_text(MINIMAL_PDB)
    return input_path


def test_repair_structure_dock_mode_adds_missing_atoms_without_hydrogens(tmp_path):
    input_path = _write_input(tmp_path)
    output_path = tmp_path / "output_dock.pdb"

    repair_structure(input_path, output_path, ph=None)

    text = output_path.read_text()
    assert " CB " in text  # ALAの側鎖が補完されている
    assert not any(line[76:78].strip() == "H" for line in text.splitlines() if line.startswith("ATOM"))


def test_repair_structure_md_mode_adds_hydrogens(tmp_path):
    input_path = _write_input(tmp_path)
    output_path = tmp_path / "output_md.pdb"

    repair_structure(input_path, output_path, ph=7.0)

    text = output_path.read_text()
    assert any(line[76:78].strip() == "H" for line in text.splitlines() if line.startswith("ATOM"))

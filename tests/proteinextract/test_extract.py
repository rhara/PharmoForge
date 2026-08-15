from structio import parse_structure

from proteinextract import extract_structure


def _atom_line(serial: int, resseq: int, chain: str, resname: str, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:>5} {'CA':<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {'C':>2}"
    )


def _water_line(serial: int, resseq: int, chain: str, x: float, y: float, z: float) -> str:
    return (
        f"HETATM{serial:>5} {'O':<4} {'HOH':>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {'O':>2}"
    )


def _write_pdb(path):
    lines = [
        _atom_line(1, 1, "A", "GLY", 0.0, 0.0, 0.0),
        _atom_line(2, 2, "A", "GLY", 1.0, 0.0, 0.0),
        _atom_line(3, 1, "B", "GLY", 0.0, 1.0, 0.0),
        _water_line(4, 100, "A", 5.0, 5.0, 5.0),
        _water_line(5, 101, "B", 6.0, 6.0, 6.0),
        "TER",
        "END",
    ]
    path.write_text("\n".join(lines) + "\n")


def test_extract_structure_filters_chains(tmp_path):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    _write_pdb(input_path)

    extract_structure(input_path, output_path, chains=["A"])

    result = parse_structure(output_path)
    assert set(result.getChids()) == {"A"}
    assert result.numAtoms() == 3  # chain Aの原子(GLY x2 + HOH x1)


def test_extract_structure_removes_water(tmp_path):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    _write_pdb(input_path)

    extract_structure(input_path, output_path, remove_water=True)

    result = parse_structure(output_path)
    assert result.numAtoms() == 3  # GLY x2 (chain A) + GLY x1 (chain B)
    water = result.select("water")
    assert water is None


def test_extract_structure_filters_chains_and_removes_water(tmp_path):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    _write_pdb(input_path)

    extract_structure(input_path, output_path, chains=["A"], remove_water=True)

    result = parse_structure(output_path)
    assert set(result.getChids()) == {"A"}
    assert result.numAtoms() == 2  # chain AのGLY x2のみ


def test_extract_structure_no_filter_keeps_everything(tmp_path):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    _write_pdb(input_path)

    extract_structure(input_path, output_path)

    result = parse_structure(output_path)
    assert result.numAtoms() == 5


def test_extract_structure_raises_when_selection_matches_nothing(tmp_path):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    _write_pdb(input_path)

    try:
        extract_structure(input_path, output_path, chains=["Z"])
        assert False, "ValueErrorが発生しなかった"
    except ValueError:
        pass


def test_extract_structure_converts_format(tmp_path):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.cif"
    _write_pdb(input_path)

    extract_structure(input_path, output_path, chains=["A"])

    assert output_path.exists()
    result = parse_structure(output_path)
    assert result.numAtoms() == 3

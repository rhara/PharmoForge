from structio import parse_structure, write_structure


def _atom_line(serial: int, resseq: int, chain: str, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:>5} {'CA':<4} {'GLY':>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {'C':>2}"
    )


def _write_pdb(path):
    lines = [_atom_line(1, 1, "A", 0.0, 0.0, 0.0), _atom_line(2, 2, "A", 1.0, 0.0, 0.0), "TER", "END"]
    path.write_text("\n".join(lines) + "\n")


def test_parse_structure_dispatches_on_extension_pdb(tmp_path):
    path = tmp_path / "input.pdb"
    _write_pdb(path)

    structure = parse_structure(path)

    assert structure.numAtoms() == 2


def test_write_structure_pdb_then_cif_roundtrip(tmp_path):
    input_path = tmp_path / "input.pdb"
    _write_pdb(input_path)
    structure = parse_structure(input_path)

    pdb_out = tmp_path / "out.pdb"
    cif_out = tmp_path / "out.cif"
    write_structure(structure, pdb_out)
    write_structure(structure, cif_out)

    assert pdb_out.exists()
    assert cif_out.exists()
    assert parse_structure(pdb_out).numAtoms() == 2
    assert parse_structure(cif_out).numAtoms() == 2


def test_write_structure_creates_parent_directories(tmp_path):
    input_path = tmp_path / "input.pdb"
    _write_pdb(input_path)
    structure = parse_structure(input_path)

    output_path = tmp_path / "nested" / "dir" / "out.pdb"
    write_structure(structure, output_path)

    assert output_path.exists()

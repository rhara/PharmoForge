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


def _write_cif_with_distinct_label_and_auth_chains(path):
    # label_asym_id(A, B)がauth_asym_id(いずれもA)と異なる最小限のmmCIF
    # (実際のPDB depositionで水分子・リガンド等が別labelチェーンに分かれるケースを模す)。
    path.write_text(
        "data_test\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM   1 C CA . GLY A 1 1 ? 0.000 0.000 0.000 1.00 40.00 ? 1 GLY A CA 1\n"
        "ATOM   2 C CA . GLY B 2 1 ? 1.000 0.000 0.000 1.00 40.00 ? 2 GLY A CA 1\n"
        "#\n"
    )


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


def test_parse_structure_cif_uses_auth_asym_id_as_chain(tmp_path):
    path = tmp_path / "input.cif"
    _write_cif_with_distinct_label_and_auth_chains(path)

    structure = parse_structure(path)

    # label_asym_idはA/Bと異なるが、auth_asym_idは両原子ともAのため1チェーンにまとまる
    assert structure.numAtoms() == 2
    assert sorted(set(structure.getChids())) == ["A"]


def test_write_structure_creates_parent_directories(tmp_path):
    input_path = tmp_path / "input.pdb"
    _write_pdb(input_path)
    structure = parse_structure(input_path)

    output_path = tmp_path / "nested" / "dir" / "out.pdb"
    write_structure(structure, output_path)

    assert output_path.exists()

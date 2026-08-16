from structio import parse_structure
from seqextract import get_chain_sequences


def _atom_line(serial, resname, resseq, chain, atom_name="CA"):
    return (
        f"ATOM  {serial:>5} {atom_name:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {atom_name[0]:>2}"
    )


def _hetatm_water_line(serial, resseq, chain):
    return (
        f"HETATM{serial:>5}  O   HOH {chain}{resseq:>4}    "
        f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          O"
    )


def _write_pdb(path):
    lines = [
        _atom_line(1, "MET", 1, "A"),
        _atom_line(2, "GLY", 2, "A"),
        _atom_line(3, "ALA", 3, "A"),
        "TER",
        _hetatm_water_line(4, 1, "B"),
        "TER",
        "END",
    ]
    path.write_text("\n".join(lines) + "\n")


def test_get_chain_sequences_extracts_protein_chain_only(tmp_path):
    path = tmp_path / "input.pdb"
    _write_pdb(path)

    chains = get_chain_sequences(parse_structure(path))

    assert [c.chain_id for c in chains] == ["A"]
    assert chains[0].sequence == "MGA"
    assert chains[0].resnums == [1, 2, 3]
    assert chains[0].length == 3


def test_get_chain_sequences_skips_non_protein_chain(tmp_path):
    path = tmp_path / "water_only.pdb"
    path.write_text(_hetatm_water_line(1, 1, "B") + "\nEND\n")

    chains = get_chain_sequences(parse_structure(path))

    assert chains == []

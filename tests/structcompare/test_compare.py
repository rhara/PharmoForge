from structio import parse_structure
from structcompare import ResidueSubstitution, match_chains, find_substitutions

_COMMON_RESNAMES = ["GLY", "GLY", "GLY", "GLY", "ALA", "GLY", "GLY", "GLY", "GLY", "GLY"]


def _atom_line(serial, resname, resseq, chain, x=0.0):
    return (
        f"ATOM  {serial:>5} {'CA':<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {'C':>2}"
    )


def _write_pdb(path, resnames, start=1, chain="A"):
    lines = [_atom_line(i + 1, rn, start + i, chain, x=float(i)) for i, rn in enumerate(resnames)]
    lines += ["TER", "END"]
    path.write_text("\n".join(lines) + "\n")


def test_match_chains_reports_identity_and_overlap(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mut_path = tmp_path / "mut.pdb"
    mutated = list(_COMMON_RESNAMES)
    mutated[4] = "SER"  # 残基番号5をALA->SERに置換
    _write_pdb(ref_path, _COMMON_RESNAMES)
    _write_pdb(mut_path, mutated)

    matches = match_chains(parse_structure(ref_path), parse_structure(mut_path))

    assert len(matches) == 1
    assert matches[0].chain_id_a == "A"
    assert matches[0].chain_id_b == "A"
    assert matches[0].seqid == 90.0
    assert matches[0].overlap == 100.0


def test_find_substitutions_detects_point_mutation_by_resnum(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mut_path = tmp_path / "mut.pdb"
    mutated = list(_COMMON_RESNAMES)
    mutated[4] = "SER"
    _write_pdb(ref_path, _COMMON_RESNAMES)
    _write_pdb(mut_path, mutated)

    report = find_substitutions(parse_structure(ref_path), parse_structure(mut_path))

    assert report.matched is True
    assert report.chain_id_b == "A"
    assert report.seqid == 90.0
    assert report.substitutions == [ResidueSubstitution(resnum=5, resname_a="ALA", resname_b="SER")]


def test_find_substitutions_empty_for_identical_sequences(tmp_path):
    p1 = tmp_path / "a.pdb"
    p2 = tmp_path / "b.pdb"
    _write_pdb(p1, _COMMON_RESNAMES)
    _write_pdb(p2, _COMMON_RESNAMES)

    report = find_substitutions(parse_structure(p1), parse_structure(p2))

    assert report.matched is True
    assert report.substitutions == []


def test_find_substitutions_not_matched_when_resnums_disjoint(tmp_path):
    p1 = tmp_path / "a.pdb"
    p2 = tmp_path / "b.pdb"
    _write_pdb(p1, _COMMON_RESNAMES, start=1)
    _write_pdb(p2, _COMMON_RESNAMES, start=100)

    report = find_substitutions(parse_structure(p1), parse_structure(p2))

    assert report.matched is False
    assert report.chain_id_b is None
    assert report.seqid is None
    assert report.substitutions == []

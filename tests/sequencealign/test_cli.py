from click.testing import CliRunner

from sequencealign.cli import sequence_align_cmd

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


def test_sequence_align_prints_report_to_stdout(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p2 = tmp_path / "mut.pdb"
    _write_pdb(p1, _COMMON_RESNAMES)
    _write_pdb(p2, _COMMON_RESNAMES)

    result = CliRunner().invoke(sequence_align_cmd, [str(p1), str(p2)])

    assert result.exit_code == 0, result.output
    assert "== Sequences (FASTA, observed residues only) ==" in result.output
    assert "== Pairwise identity ==" in result.output
    assert "== Substitutions relative to reference ==" not in result.output


def test_sequence_align_with_reference_structure(tmp_path):
    p1 = tmp_path / "ref.pdb"
    p2 = tmp_path / "mut.pdb"
    _write_pdb(p1, _COMMON_RESNAMES)
    _write_pdb(p2, _COMMON_RESNAMES)

    result = CliRunner().invoke(sequence_align_cmd, [str(p1), str(p2), "--reference", "ref:A"])

    assert result.exit_code == 0, result.output
    assert "== Substitutions relative to reference ==" in result.output


def test_sequence_align_indir_resolves_stems(tmp_path):
    (tmp_path / "a.pdb").write_text("\n".join([_atom_line(1, "GLY", 1, "A"), "TER", "END"]) + "\n")

    result = CliRunner().invoke(sequence_align_cmd, ["--indir", str(tmp_path), "a"])

    assert result.exit_code == 0, result.output
    assert ">a:A" in result.output


def test_sequence_align_writes_output_file(tmp_path):
    p1 = tmp_path / "ref.pdb"
    _write_pdb(p1, _COMMON_RESNAMES)
    out_path = tmp_path / "report.txt"

    result = CliRunner().invoke(sequence_align_cmd, [str(p1), "-o", str(out_path)])

    assert result.exit_code == 0, result.output
    assert out_path.exists()
    assert "== Sequences (FASTA, observed residues only) ==" in out_path.read_text()


def test_sequence_align_requires_at_least_one_path():
    result = CliRunner().invoke(sequence_align_cmd, [])
    assert result.exit_code != 0

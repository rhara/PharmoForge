import pytest

from sequencealign.report import (
    build_report,
    format_fasta,
    format_identity_matrix,
    format_mutation_report,
    load_labeled_structures,
)

_COMMON_RESNAMES = ["GLY", "GLY", "GLY", "GLY", "ALA", "GLY", "GLY", "GLY", "GLY", "GLY"]
_COMMON_SEQUENCE = "GGGGAGGGGG"

# アラインメントの曖昧さを避けるため、ギャップ検出のテストには非反復配列を使う(MENFQKVEKI)
_VARIED_RESNAMES = ["MET", "GLU", "ASN", "PHE", "GLN", "LYS", "VAL", "GLU", "LYS", "ILE"]
_VARIED_SEQUENCE = "MENFQKVEKI"


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


@pytest.fixture
def two_structures(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mut_path = tmp_path / "mut.pdb"
    mutated = list(_COMMON_RESNAMES)
    mutated[4] = "SER"  # 残基番号5をALA->SERに置換
    _write_pdb(ref_path, _COMMON_RESNAMES)
    _write_pdb(mut_path, mutated)
    return load_labeled_structures([ref_path, mut_path])


@pytest.fixture
def structures_with_gap(tmp_path):
    # refは残基1-10(MENFQKVEKI)、gappedは残基6-9(KVEK)のみ保持(1-5, 10が欠損)
    ref_path = tmp_path / "ref.pdb"
    gapped_path = tmp_path / "gapped.pdb"
    _write_pdb(ref_path, _VARIED_RESNAMES, start=1)
    _write_pdb(gapped_path, _VARIED_RESNAMES[5:9], start=6)
    return load_labeled_structures([ref_path, gapped_path])


def test_format_fasta_lists_all_chains(two_structures):
    fasta = format_fasta(two_structures)

    assert ">ref:A length=10 range=1-10" in fasta
    assert ">mut:A length=10 range=1-10" in fasta
    assert _COMMON_SEQUENCE in fasta


def test_format_identity_matrix_reports_pairwise_identity(two_structures):
    matrix = format_identity_matrix(two_structures)

    assert "ref:A  vs  mut:A" in matrix
    assert "identity= 90.0%" in matrix


def test_format_mutation_report_vs_structure_reference(two_structures):
    report = format_mutation_report(two_structures, "ref:A")

    assert "基準: ref:A" in report
    assert "mut: 1箇所" in report
    assert "A5S" in report


def test_format_mutation_report_vs_sequence_reference(two_structures):
    report = format_mutation_report(two_structures, _COMMON_SEQUENCE)

    assert "基準配列: 10残基" in report
    assert "ref:A: 置換なし" in report
    assert "mut:A: 1箇所" in report
    assert "A5S" in report
    assert "構造残基番号=5" in report


def test_format_mutation_report_rejects_invalid_sequence(two_structures):
    with pytest.raises(ValueError):
        format_mutation_report(two_structures, "not-a-sequence-or-ref")


def test_format_mutation_report_unknown_structure_label_raises(two_structures):
    with pytest.raises(ValueError):
        format_mutation_report(two_structures, "nonexistent:A")


def test_build_report_includes_all_sections(two_structures):
    report = build_report(two_structures, reference="ref:A")

    assert "== 配列(FASTA、観測された残基のみ) ==" in report
    assert "== Pairwise identity ==" in report
    assert "== 基準配列に対する置換 ==" in report


def test_build_report_omits_mutation_section_without_reference(two_structures):
    report = build_report(two_structures, reference=None)

    assert "== 基準配列に対する置換 ==" not in report


def test_format_mutation_report_vs_structure_shows_gap(structures_with_gap):
    report = format_mutation_report(structures_with_gap, "ref:A")

    assert "gapped: 置換なし" in report
    assert "欠損: 基準のみ(対象で欠損): 1-5, 10" in report


def test_format_mutation_report_vs_sequence_shows_gap(structures_with_gap):
    report = format_mutation_report(structures_with_gap, _VARIED_SEQUENCE)

    assert "gapped:A: 置換なし" in report
    assert "欠損: 基準1-5(5残基), 基準10(1残基)" in report

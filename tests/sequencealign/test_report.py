import pytest

from sequencealign.report import (
    _format_ruler,
    build_report,
    format_alignment_block,
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

    assert "reference: ref:A" in report
    assert "mut: 1 substitution(s)" in report
    assert "A5S" in report


def test_format_mutation_report_vs_sequence_reference(two_structures):
    report = format_mutation_report(two_structures, _COMMON_SEQUENCE)

    assert "reference sequence: 10 residues" in report
    assert "ref:A: no substitutions" in report
    assert "mut:A: 1 substitution(s)" in report
    assert "A5S" in report
    assert "structure resnum=5" in report


def test_format_mutation_report_rejects_invalid_sequence(two_structures):
    with pytest.raises(ValueError):
        format_mutation_report(two_structures, "not-a-sequence-or-ref")


def test_format_mutation_report_unknown_structure_label_raises(two_structures):
    with pytest.raises(ValueError):
        format_mutation_report(two_structures, "nonexistent:A")


def test_build_report_includes_all_sections(two_structures):
    report = build_report(two_structures, reference="ref:A")

    assert "== Sequences (FASTA, observed residues only) ==" in report
    assert "== Pairwise identity ==" in report
    assert "== Substitutions relative to reference ==" in report


def test_build_report_omits_mutation_section_without_reference(two_structures):
    report = build_report(two_structures, reference=None)

    assert "== Substitutions relative to reference ==" not in report


def test_format_mutation_report_vs_structure_shows_gap(structures_with_gap):
    report = format_mutation_report(structures_with_gap, "ref:A")

    assert "gapped: no substitutions" in report
    assert "gaps: reference only (missing in target): 1-5, 10" in report


def test_format_mutation_report_vs_sequence_shows_gap(structures_with_gap):
    report = format_mutation_report(structures_with_gap, _VARIED_SEQUENCE)

    assert "gapped:A: no substitutions" in report
    assert "gaps: ref 1-5 (5 residue(s)), ref 10 (1 residue(s))" in report


def test_format_alignment_block_lines_up_by_resnum(structures_with_gap):
    block = format_alignment_block(structures_with_gap)

    assert "-- 1-10 --" in block
    assert "ref:A     MENFQKVEKI" in block
    assert "gapped:A  -----KVEK-" in block


def test_format_ruler_places_number_and_tick_at_multiples_of_ten():
    numbers, ticks = _format_ruler(101, 100)

    assert numbers[7:10] == "110"  # resnum 110 は列9(0始まり)、末尾の'0'がその列に来る
    assert ticks[9] == "|"
    assert numbers[97:100] == "200"  # resnum 200 は列99(0始まり)
    assert ticks[99] == "|"
    assert ticks.replace("|", "").strip() == ""


def test_format_ruler_no_tick_for_short_block():
    numbers, ticks = _format_ruler(1, 5)

    assert ticks == " " * 5
    assert numbers == " " * 5


def test_format_ruler_omits_tick_that_would_be_truncated_at_left_edge():
    # start_resnum=449: 最初の10の倍数(450)は列1で、"450"の上位桁がブロック外にはみ出す。
    # 上位桁が欠けて'0'だけが見える(実際の値を誤読させる)くらいなら目盛りごと省略する。
    numbers, ticks = _format_ruler(449, 100)

    assert "450" not in numbers
    assert numbers[:9] == " " * 9  # resnum 450(列1)は省略され、次の460(列11)から始まる
    assert numbers[9:12] == "460"
    assert ticks[11] == "|"
    assert ticks[:11] == " " * 11


def test_format_alignment_block_includes_ruler(structures_with_gap):
    block = format_alignment_block(structures_with_gap)

    assert "          10" in block  # 位置10の目盛り(数字行)
    assert "           |" in block  # 位置10の目盛り(縦棒行)


def test_format_alignment_block_wraps_at_given_width():
    # 25残基を10残基/行で折り返す
    resnames = ["GLY"] * 25
    import tempfile
    from pathlib import Path

    d = Path(tempfile.mkdtemp())
    path = d / "a.pdb"
    _write_pdb(path, resnames)
    structures = load_labeled_structures([path])

    block = format_alignment_block(structures, width=10)

    assert "-- 1-10 --" in block
    assert "-- 11-20 --" in block
    assert "-- 21-25 --" in block


def test_format_alignment_block_empty_when_no_chains(tmp_path):
    path = tmp_path / "water.pdb"
    path.write_text(
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nEND\n"
    )
    structures = load_labeled_structures([path])

    assert format_alignment_block(structures) == "(no protein chains found)\n"


def test_build_report_includes_alignment_block(two_structures):
    report = build_report(two_structures, reference=None)

    assert "== Alignment (by residue number) ==" in report

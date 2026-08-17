import pytest

from sequencealign.report import (
    _format_ruler,
    build_report,
    format_alignment_block,
    format_alignment_block_by_sequence,
    format_coverage_matrix,
    format_fasta,
    format_identity_coverage_matrix,
    format_identity_matrix,
    format_mutation_report,
    load_labeled_structures,
)

_COMMON_RESNAMES = ["GLY", "GLY", "GLY", "GLY", "ALA", "GLY", "GLY", "GLY", "GLY", "GLY"]
_COMMON_SEQUENCE = "GGGGAGGGGG"

# アラインメントの曖昧さを避けるため、ギャップ検出のテストには非反復配列を使う(MENFQKVEKI)
_VARIED_RESNAMES = ["MET", "GLU", "ASN", "PHE", "GLN", "LYS", "VAL", "GLU", "LYS", "ILE"]
_VARIED_SEQUENCE = "MENFQKVEKI"

# P0DTD1(全長polyprotein)/6LU7(ローカル採番のMproドメインのみ)を模した組み合わせ:
# 全長配列の中に非反復ドメイン配列を埋め込み、ドメイン側は独立した1始まりの残基番号を持つ。
_DOMAIN_PREFIX = "AAAAA"
_DOMAIN_SUFFIX = "CCCCC"
_FULL_LENGTH_SEQUENCE = _DOMAIN_PREFIX + _VARIED_SEQUENCE + _DOMAIN_SUFFIX  # 20残基、ドメインは6-15


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


@pytest.fixture
def full_length_and_domain_only(tmp_path):
    # P0DTD1(全長polyprotein FASTA)/6LU7(ローカル採番のMproドメインのみのPDB)を模す:
    # full_pdbは全長配列(1-20)、domain_pdbは同じドメイン配列だが独立採番(1-10)を持つ。
    fasta_path = tmp_path / "full.fasta"
    fasta_path.write_text(f">full\n{_FULL_LENGTH_SEQUENCE}\n")
    domain_path = tmp_path / "domain.pdb"
    _write_pdb(domain_path, _VARIED_RESNAMES, start=1)
    return load_labeled_structures([fasta_path, domain_path])


def test_format_fasta_lists_all_chains(two_structures):
    fasta = format_fasta(two_structures)

    assert ">ref:A length=10 range=1-10" in fasta
    assert ">mut:A length=10 range=1-10" in fasta
    assert _COMMON_SEQUENCE in fasta


def test_format_identity_matrix_reports_pairwise_identity(two_structures):
    matrix = format_identity_matrix(two_structures)

    lines = matrix.splitlines()
    assert lines[0].split() == ["ref:A", "mut:A"]  # ヘッダー行(列ラベル)
    ref_row, mut_row = lines[1].split(), lines[2].split()
    assert ref_row == ["ref:A", "-", "90.0"]
    assert mut_row == ["mut:A", "90.0", "-"]  # 対称なグリッド


def test_format_mutation_report_vs_structure_reference(two_structures):
    report = format_mutation_report(two_structures, "ref:A")

    assert "reference: ref:A" in report
    assert _COMMON_SEQUENCE in report  # 基準チェーンの配列自体も出力される
    assert "mut: 1 substitution(s)" in report
    assert "A5S" in report


def test_format_mutation_report_vs_sequence_reference(two_structures):
    report = format_mutation_report(two_structures, _COMMON_SEQUENCE)

    assert "reference sequence: 10 residues" in report
    assert _COMMON_SEQUENCE in report  # 指定した基準配列自体も出力される
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


def test_build_report_includes_pairwise_identity_and_alignment_only(two_structures):
    report = build_report(two_structures)

    assert "== Pairwise identity/coverage ==" in report  # identity_format既定は"combined"
    assert "== Alignment (sequence-aligned) ==" in report  # method既定は"align"
    # FASTA・Substitutionsセクションはユーザー要望により出力しない
    assert "== Sequences" not in report
    assert "== Substitutions" not in report


def test_build_report_includes_fasta_input_as_alignment_row(two_structures, tmp_path):
    fasta_path = tmp_path / "reference.fasta"
    fasta_path.write_text(f">reference\n{_COMMON_SEQUENCE}\n")
    structures = two_structures + load_labeled_structures([fasta_path])

    report = build_report(structures)

    assert "reference:A" in report


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

    assert "ref:A     MENFQKVEKI" in block
    assert "gapped:A  -----KVEK-" in block


def test_format_alignment_block_includes_fasta_input_row(structures_with_gap, tmp_path):
    fasta_path = tmp_path / "reference.fasta"
    fasta_path.write_text(f">reference\n{_VARIED_SEQUENCE}\n")
    structures = structures_with_gap + load_labeled_structures([fasta_path])

    block = format_alignment_block(structures)

    assert "reference:A  MENFQKVEKI" in block
    assert "ref:A        MENFQKVEKI" in block
    assert "gapped:A     -----KVEK-" in block


def test_format_alignment_block_misaligns_mismatched_numbering(full_length_and_domain_only):
    # 番号ベースでは、ドメイン側のローカル採番(1-10)が全長配列の先頭(1-10、AAAAAMENFQ)と
    # そのまま重なってしまい、ドメイン本来の位置(6-15)には並ばない(既知の制約)。
    block = format_alignment_block(full_length_and_domain_only)

    assert "full:A    AAAAAMENFQKVEKICCCCC" in block
    assert "domain:A  MENFQKVEKI----------" in block


def test_format_alignment_block_by_sequence_aligns_domain_at_correct_offset(full_length_and_domain_only):
    # 配列アラインメントベースでは、ドメインの配列が全長配列内の本来の位置(6-15)に
    # 正しく並ぶ(先頭の構造=fullが基準)。
    block = format_alignment_block_by_sequence(full_length_and_domain_only)

    assert "full:A    AAAAAMENFQKVEKICCCCC" in block
    assert "domain:A  -----MENFQKVEKI-----" in block


def test_format_alignment_block_by_sequence_empty_when_no_chains(tmp_path):
    path = tmp_path / "water.pdb"
    path.write_text(
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nEND\n"
    )
    structures = load_labeled_structures([path])

    assert format_alignment_block_by_sequence(structures) == "(no protein chains found)\n"


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


def test_format_ruler_left_aligns_number_that_would_overflow_left_edge():
    # start_resnum=449: 最初の10の倍数(450)は列1で、右揃えだと"450"の上位桁が
    # ブロック外にはみ出す。上位桁が欠けて'0'だけが見える(実際の値を誤読させる)ことを
    # 避けるため、列0に寄せて数字全体("450")を表示する。tick自体は本来の列(1)のまま。
    numbers, ticks = _format_ruler(449, 100)

    assert numbers[0:3] == "450"
    assert ticks[1] == "|"
    assert numbers[9:12] == "460"
    assert ticks[11] == "|"


def test_load_labeled_structures_reads_single_record_fasta(tmp_path):
    fasta_path = tmp_path / "P0DTD1.fasta"
    fasta_path.write_text(">sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab\nMENFQKVEKI\n")

    structures = load_labeled_structures([fasta_path])

    assert len(structures) == 1
    s = structures[0]
    assert s.label == "P0DTD1"
    assert s.atoms is None
    assert len(s.chains) == 1
    assert s.chains[0].chain_id == "A"
    assert s.chains[0].sequence == "MENFQKVEKI"
    assert s.chains[0].resnums == list(range(1, 11))


def test_load_labeled_structures_reads_multi_record_fasta(tmp_path):
    fasta_path = tmp_path / "complex.fasta"
    fasta_path.write_text(">chain1\nGGGGAGGGGG\n>chain2\nMENFQKVEKI\n")

    structures = load_labeled_structures([fasta_path])

    assert len(structures) == 1
    chains = structures[0].chains
    assert [c.chain_id for c in chains] == ["A", "B"]
    assert chains[0].sequence == "GGGGAGGGGG"
    assert chains[1].sequence == "MENFQKVEKI"


def test_format_identity_matrix_includes_fasta_only_chains(two_structures, tmp_path):
    fasta_path = tmp_path / "seqonly.fasta"
    fasta_path.write_text(f">seqonly\n{_COMMON_SEQUENCE}\n")
    structures = two_structures + load_labeled_structures([fasta_path])

    matrix = format_identity_matrix(structures)

    # atomsを持たないFASTA由来のチェーンも、配列アラインメントベースでグリッドに含まれる。
    assert "seqonly:A" in matrix
    lines = matrix.splitlines()
    assert lines[0].split() == ["ref:A", "mut:A", "seqonly:A"]
    ref_row = lines[1].split()
    assert ref_row == ["ref:A", "-", "90.0", "100.0"]  # ref=_COMMON_SEQUENCEに1残基置換したものがmut


def test_format_identity_matrix_fewer_than_two_chains(tmp_path):
    path = tmp_path / "single.fasta"
    path.write_text(">a\nMENFQKVEKI\n")
    structures = load_labeled_structures([path])

    assert format_identity_matrix(structures) == "(fewer than two protein chains found)\n"


def test_format_coverage_matrix_is_asymmetric(full_length_and_domain_only):
    # full(20残基)を基準にすると、domain(10残基)はその半分しかカバーしない(50.0%)。
    # domain(10残基)を基準にすると、full内の対応部分は全て含まれる(100.0%)。
    matrix = format_coverage_matrix(full_length_and_domain_only)

    lines = matrix.splitlines()
    assert lines[0].split() == ["full:A", "domain:A"]
    full_row = lines[1].split()
    domain_row = lines[2].split()
    assert full_row == ["full:A", "-", "50.0"]
    assert domain_row == ["domain:A", "100.0", "-"]


def test_format_coverage_matrix_fewer_than_two_chains(tmp_path):
    path = tmp_path / "single.fasta"
    path.write_text(">a\nMENFQKVEKI\n")
    structures = load_labeled_structures([path])

    assert format_coverage_matrix(structures) == "(fewer than two protein chains found)\n"


def test_format_identity_coverage_matrix_combines_both_values(full_length_and_domain_only):
    matrix = format_identity_coverage_matrix(full_length_and_domain_only)

    lines = matrix.splitlines()
    full_row = lines[1].split()
    domain_row = lines[2].split()
    assert full_row == ["full:A", "-", "100.0/50.0"]
    assert domain_row == ["domain:A", "100.0/100.0", "-"]


def test_build_report_identity_format_separate_shows_two_tables(full_length_and_domain_only):
    report = build_report(full_length_and_domain_only, identity_format="separate")

    assert "== Pairwise identity ==" in report
    assert "== Coverage ==" in report
    assert "== Pairwise identity/coverage ==" not in report
    assert "100.0/50.0" not in report  # combined形式のセル表記は出ない


def test_build_report_identity_format_combined_shows_one_table(full_length_and_domain_only):
    report = build_report(full_length_and_domain_only, identity_format="combined")

    assert "== Pairwise identity/coverage ==" in report
    assert "== Coverage ==" not in report
    assert "100.0/50.0" in report


def test_format_mutation_report_vs_structure_rejects_fasta_reference(two_structures, tmp_path):
    fasta_path = tmp_path / "seqonly.fasta"
    fasta_path.write_text(f">seqonly\n{_COMMON_SEQUENCE}\n")
    structures = two_structures + load_labeled_structures([fasta_path])

    with pytest.raises(ValueError):
        format_mutation_report(structures, "seqonly:A")


def test_format_mutation_report_vs_structure_notes_fasta_only_target(two_structures, tmp_path):
    fasta_path = tmp_path / "seqonly.fasta"
    fasta_path.write_text(f">seqonly\n{_COMMON_SEQUENCE}\n")
    structures = two_structures + load_labeled_structures([fasta_path])

    report = format_mutation_report(structures, "ref:A")

    assert "seqonly: no atomic structure (loaded from FASTA)" in report


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

    # 25残基が10残基/行で3ブロック(空行区切り)に折り返される
    assert len(block.strip("\n").split("\n\n")) == 3
    assert "a:A  GGGGGGGGGG" in block
    assert "a:A  GGGGG\n" in block


def test_format_alignment_block_empty_when_no_chains(tmp_path):
    path = tmp_path / "water.pdb"
    path.write_text(
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nEND\n"
    )
    structures = load_labeled_structures([path])

    assert format_alignment_block(structures) == "(no protein chains found)\n"


def test_build_report_includes_alignment_block(two_structures):
    report = build_report(two_structures)

    assert "== Alignment (sequence-aligned) ==" in report  # method既定は"align"


def test_build_report_method_number_uses_residue_number_alignment(two_structures):
    report = build_report(two_structures, method="number")

    assert "== Alignment (by residue number) ==" in report


def test_build_report_method_align_uses_sequence_alignment(full_length_and_domain_only):
    report = build_report(full_length_and_domain_only, method="align")

    assert "== Alignment (sequence-aligned) ==" in report
    assert "domain:A  -----MENFQKVEKI-----" in report

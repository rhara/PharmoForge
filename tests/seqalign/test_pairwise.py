from seqalign import align_to_reference

_REF = "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNH"


def test_align_to_reference_finds_identical_sequence():
    result = align_to_reference(_REF, _REF, list(range(1, len(_REF) + 1)))
    assert result.substitutions == []
    assert result.identity == 100.0
    assert result.coverage == 100.0
    assert result.aligned_length == len(_REF)
    assert result.query_by_ref_pos == {i + 1: c for i, c in enumerate(_REF)}


def test_align_to_reference_query_by_ref_pos_places_domain_at_correct_offset():
    # queryは基準配列の一部(11-20番目、1始まり)のみに対応する短いドメイン配列とみなす
    domain = _REF[10:20]
    result = align_to_reference(_REF, domain, list(range(1, len(domain) + 1)))

    assert result.query_by_ref_pos == {i + 11: c for i, c in enumerate(domain)}
    assert 1 not in result.query_by_ref_pos
    assert 21 not in result.query_by_ref_pos


def test_align_to_reference_finds_point_substitution():
    query = "X" + _REF[1:]  # 先頭を置換
    result = align_to_reference(_REF, query, list(range(1, len(query) + 1)))

    assert len(result.substitutions) == 1
    sub = result.substitutions[0]
    assert sub.ref_pos == 1
    assert sub.ref_aa == "M"
    assert sub.query_aa == "X"
    assert sub.query_resnum == 1
    assert result.identity == 100.0 * (len(_REF) - 1) / len(_REF)


def test_align_to_reference_maps_resnums_across_gap():
    # queryは先頭15残基が欠損(構造の残基番号は16始まり)しているとみなす
    query = _REF[15:]
    resnums = list(range(16, 16 + len(query)))

    result = align_to_reference(_REF, query, resnums)

    assert result.substitutions == []
    assert result.identity == 100.0
    assert result.aligned_length == len(query)
    assert result.coverage == 100.0 * len(query) / len(_REF)
    assert len(result.gaps) == 1
    gap = result.gaps[0]
    assert gap.kind == "deletion"
    assert gap.ref_start == 1
    assert gap.ref_end == 15
    assert gap.length == 15
    assert gap.before_query_resnum is None  # N末端側の欠損
    assert gap.after_query_resnum == 16


def test_align_to_reference_detects_internal_deletion():
    # 残基21-30(1始まり)が構造上欠損しているとみなす
    query = _REF[:20] + _REF[30:]
    resnums = list(range(1, 21)) + list(range(31, 31 + len(query) - 20))

    result = align_to_reference(_REF, query, resnums)

    assert result.substitutions == []
    deletions = [g for g in result.gaps if g.kind == "deletion"]
    assert len(deletions) == 1
    gap = deletions[0]
    assert gap.ref_start == 21
    assert gap.ref_end == 30
    assert gap.length == 10
    assert gap.before_query_resnum == 20
    assert gap.after_query_resnum == 31


def test_align_to_reference_detects_insertion():
    # queryの先頭に基準配列にない5残基(発現タグ等)が付加されているとみなす
    tag = "AAAAA"
    query = tag + _REF
    resnums = list(range(-4, 1)) + list(range(1, len(_REF) + 1))

    result = align_to_reference(_REF, query, resnums)

    assert result.substitutions == []
    insertions = [g for g in result.gaps if g.kind == "insertion"]
    assert len(insertions) == 1
    gap = insertions[0]
    assert gap.length == 5
    assert gap.ref_start is None
    assert gap.before_query_resnum is None  # N末端側なので直前の文脈はない
    assert gap.after_query_resnum == 1


def test_align_to_reference_maps_resnum_correctly_after_gap_with_substitution():
    # 残基11-15(1始まり)が構造上欠損しているとみなし、離れた位置(残基40)を置換
    query_chars = list(_REF[:10] + _REF[15:])
    query_chars[34] = "X"  # 元のREF[39] (1始まりで残基40) に対応する位置
    query = "".join(query_chars)
    resnums = list(range(1, 11)) + list(range(16, 16 + len(query) - 10))

    result = align_to_reference(_REF, query, resnums)

    assert len(result.substitutions) == 1
    sub = result.substitutions[0]
    assert sub.ref_pos == 40
    assert sub.query_resnum == 40
    assert sub.query_aa == "X"

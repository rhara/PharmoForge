"""任意のアミノ酸配列(1文字表記)同士のペアワイズグローバルアラインメント。

`structcompare`はProDyの構造(Atomic)同士の比較に特化しているが、構造を伴わない
任意配列(例: UniProt正規配列やユーザー指定の基準配列)と構造由来配列を比較する
場合は比較対象の一方にAtomicが存在しないため、Biopythonの
`Bio.Align.PairwiseAligner`を直接用いる。
"""

from dataclasses import dataclass

from Bio.Align import PairwiseAligner, substitution_matrices

from core.logging_utils import get_logger

logger = get_logger(__name__)

_SUBSTITUTION_MATRIX = substitution_matrices.load("BLOSUM62")
_OPEN_GAP_SCORE = -11.0
_EXTEND_GAP_SCORE = -1.0


@dataclass
class SequenceSubstitution:
    """対応する位置でアミノ酸が異なる箇所。"""

    ref_pos: int  # 基準配列内の位置(1始まり)
    ref_aa: str
    query_resnum: int  # 比較対象(構造)側の残基番号
    query_aa: str


@dataclass
class SequenceGap:
    """基準配列に対する欠失(query側に residue がない)、または挿入(query側にのみ residue がある)領域。

    PDB構造では電子密度が見えない領域(ループ等)が欠落することが多く、置換とは別に
    可視化する。`kind`は`"deletion"`(基準配列にはあるがqueryにない)または
    `"insertion"`(queryにあり基準配列にない)。`before_query_resnum`/`after_query_resnum`は
    いずれもこの領域自体の残基番号ではなく、query側でこの領域の直前・直後にある
    (アラインメントされた)残基番号(文脈)。該当がない場合(配列の末端)は`None`。
    """

    kind: str
    ref_start: int | None  # 基準配列内の欠失範囲(1始まり、insertionの場合None)
    ref_end: int | None
    length: int
    before_query_resnum: int | None  # 直前の(欠落していない)query側残基番号。末端の場合None
    after_query_resnum: int | None  # 直後の(欠落していない)query側残基番号。末端の場合None


@dataclass
class AlignmentResult:
    """`ref_sequence`と`query_sequence`のグローバルアラインメント結果。"""

    identity: float  # % (アラインメントされた位置のうち一致した割合)
    coverage: float  # % (ref_sequence全体のうちアラインメントされた位置の割合)
    aligned_length: int  # アラインメントされた位置数(ギャップを除く)
    substitutions: list[SequenceSubstitution]
    gaps: list[SequenceGap]


def align_to_reference(ref_sequence: str, query_sequence: str, query_resnums: list[int]) -> AlignmentResult:
    """`ref_sequence`と`query_sequence`(長さ`len(query_resnums)`)をグローバル
    アラインメントし、%identity・%coverage・アミノ酸置換・欠失/挿入領域を返す。

    置換は対応する位置(どちらにもギャップが入らない位置)のみを対象とする
    (挿入・欠失そのものは`gaps`側で扱う)。
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = _SUBSTITUTION_MATRIX
    aligner.open_gap_score = _OPEN_GAP_SCORE
    aligner.extend_gap_score = _EXTEND_GAP_SCORE
    alignment = aligner.align(ref_sequence, query_sequence)[0]

    ref_blocks, query_blocks = alignment.aligned
    ref_blocks = [(int(s), int(e)) for s, e in ref_blocks]
    query_blocks = [(int(s), int(e)) for s, e in query_blocks]
    substitutions = []
    gaps = []
    aligned_length = 0
    n_identical = 0
    prev_ref_end = 0
    prev_query_end = 0
    for (ref_start, ref_end), (q_start, q_end) in zip(ref_blocks, query_blocks):
        if ref_start > prev_ref_end:
            gaps.append(
                SequenceGap(
                    kind="deletion",
                    ref_start=prev_ref_end + 1,
                    ref_end=ref_start,
                    length=ref_start - prev_ref_end,
                    before_query_resnum=query_resnums[prev_query_end - 1] if prev_query_end > 0 else None,
                    after_query_resnum=query_resnums[q_start] if q_start < len(query_resnums) else None,
                )
            )
        if q_start > prev_query_end:
            gaps.append(
                SequenceGap(
                    kind="insertion",
                    ref_start=None,
                    ref_end=None,
                    length=q_start - prev_query_end,
                    before_query_resnum=query_resnums[prev_query_end - 1] if prev_query_end > 0 else None,
                    after_query_resnum=query_resnums[q_start] if q_start < len(query_resnums) else None,
                )
            )

        aligned_length += ref_end - ref_start
        for offset in range(ref_end - ref_start):
            ref_i = ref_start + offset
            q_i = q_start + offset
            if ref_sequence[ref_i] == query_sequence[q_i]:
                n_identical += 1
            else:
                substitutions.append(
                    SequenceSubstitution(
                        ref_pos=int(ref_i) + 1,
                        ref_aa=ref_sequence[ref_i],
                        query_resnum=query_resnums[q_i],
                        query_aa=query_sequence[q_i],
                    )
                )
        prev_ref_end = ref_end
        prev_query_end = q_end

    if len(ref_sequence) > prev_ref_end:
        gaps.append(
            SequenceGap(
                kind="deletion",
                ref_start=prev_ref_end + 1,
                ref_end=len(ref_sequence),
                length=len(ref_sequence) - prev_ref_end,
                before_query_resnum=query_resnums[prev_query_end - 1] if prev_query_end > 0 else None,
                after_query_resnum=None,
            )
        )
    if len(query_sequence) > prev_query_end:
        gaps.append(
            SequenceGap(
                kind="insertion",
                ref_start=None,
                ref_end=None,
                length=len(query_sequence) - prev_query_end,
                before_query_resnum=query_resnums[prev_query_end - 1] if prev_query_end > 0 else None,
                after_query_resnum=None,
            )
        )

    identity = 100.0 * n_identical / aligned_length if aligned_length else 0.0
    coverage = 100.0 * aligned_length / len(ref_sequence) if ref_sequence else 0.0
    logger.info(
        "align_to_reference: %d substitution(s), %d gap(s) found (identity=%.1f%%, coverage=%.1f%%)",
        len(substitutions), len(gaps), identity, coverage,
    )
    return AlignmentResult(
        identity=identity,
        coverage=coverage,
        aligned_length=aligned_length,
        substitutions=substitutions,
        gaps=gaps,
    )

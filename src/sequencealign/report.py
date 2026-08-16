"""複数構造間の蛋白配列比較レポート(FASTA・pairwise identity・基準配列に対する置換一覧)。"""

import re
from dataclasses import dataclass
from pathlib import Path

from prody.atomic.atomic import AAMAP, Atomic

from core.logging_utils import get_logger
from seqalign import align_to_reference
from seqextract import ChainSequence, get_chain_sequences
from structcompare import find_substitutions, match_chains
from structio import parse_structure

logger = get_logger(__name__)

_FASTA_WIDTH = 60
_ALIGN_WIDTH = 100
# 標準20種 + 曖昧/非標準コード(Asx/Glx/Xle/Sec/Pyl/不明)。matchChains側のAAMAPと合わせる。
_SEQUENCE_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZJUO]+$")


@dataclass
class LabeledStructure:
    """入力トークンから解決したファイル名(拡張子抜き)をラベルとする構造。"""

    label: str
    atoms: Atomic
    chains: list[ChainSequence]


def load_labeled_structures(paths: list[Path]) -> list[LabeledStructure]:
    """構造ファイルを読み込み、ファイル名(拡張子抜き)をラベルとして付与する。"""
    structures = []
    for path in paths:
        logger.info("Loading structure from %s ...", path)
        atoms = parse_structure(path)
        chains = get_chain_sequences(atoms)
        logger.info("  -> %d protein chain(s): %s", len(chains), [c.chain_id for c in chains])
        structures.append(LabeledStructure(label=path.stem, atoms=atoms, chains=chains))
    return structures


def _wrap(seq: str, width: int = _FASTA_WIDTH) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def format_fasta(structures: list[LabeledStructure]) -> str:
    """全構造・全蛋白チェーンの配列をFASTA形式で整形する。"""
    lines = []
    for s in structures:
        for c in s.chains:
            lines.append(f">{s.label}:{c.chain_id} length={c.length} range={c.resnums[0]}-{c.resnums[-1]}")
            lines.append(_wrap(c.sequence))
    if not lines:
        return "(no protein chains found)\n"
    return "\n".join(lines) + "\n"


def format_identity_matrix(structures: list[LabeledStructure]) -> str:
    """全構造の組み合わせについて、チェーン単位のpairwise identity/overlapを一覧化する。"""
    rows = []
    for i in range(len(structures)):
        for j in range(i + 1, len(structures)):
            sa, sb = structures[i], structures[j]
            if not sa.chains or not sb.chains:
                continue
            for m in match_chains(sa.atoms, sb.atoms):
                rows.append(
                    f"{sa.label}:{m.chain_id_a}  vs  {sb.label}:{m.chain_id_b}"
                    f"  identity={m.seqid:5.1f}%  overlap={m.overlap:5.1f}%  (n={m.n_matched})"
                )
    if not rows:
        return "(no comparable chain pairs found)\n"
    return "\n".join(rows) + "\n"


def format_alignment_block(structures: list[LabeledStructure], width: int = _ALIGN_WIDTH) -> str:
    """全構造・全蛋白チェーンの配列を、残基番号を共通の軸として横並びに整列表示する
    (`width`残基ごとに折り返す)。

    配列アラインメントは行わず、残基番号が一致する列に同じアミノ酸が並ぶ前提で
    並べる(構造間でPDBの残基番号が揃っている前提。`pf align-view --method number`
    と同じ前提)。観測されていない残基は`-`で埋める。異なる蛋白の構造を混在させると
    無意味な結果になる点に注意(通常は同一蛋白の複数構造を対象とする)。
    """
    entries = [(f"{s.label}:{c.chain_id}", c) for s in structures for c in s.chains]
    if not entries:
        return "(no protein chains found)\n"

    all_resnums = {r for _, c in entries for r in c.resnums}
    min_resnum, max_resnum = min(all_resnums), max(all_resnums)
    label_width = max(len(label) for label, _ in entries)

    rows = []
    for label, c in entries:
        seq_by_resnum = dict(zip(c.resnums, c.sequence))
        padded = "".join(seq_by_resnum.get(r, "-") for r in range(min_resnum, max_resnum + 1))
        rows.append((label, padded))

    total_length = max_resnum - min_resnum + 1
    blocks = []
    for block_start in range(0, total_length, width):
        block_end = min(block_start + width, total_length)
        block_first_resnum = min_resnum + block_start
        header = f"-- {block_first_resnum}-{min_resnum + block_end - 1} --"
        number_line, tick_line = _format_ruler(block_first_resnum, block_end - block_start)
        indent = " " * (label_width + 2)
        block_lines = [header, indent + number_line, indent + tick_line]
        block_lines += [f"{label.ljust(label_width)}  {seq[block_start:block_end]}" for label, seq in rows]
        blocks.append("\n".join(block_lines))
    return "\n\n".join(blocks) + "\n"


def _format_ruler(start_resnum: int, block_width: int) -> tuple[str, str]:
    """10残基ごとに残基番号とその位置を示す目盛り(2行: 数字の行、`|`の行)を作る。

    数字はその残基番号の列で右端が揃うように配置する(例: 残基120の場合、'0'が
    resnum=120の列に来る)。
    """
    numbers = [" "] * block_width
    ticks = [" "] * block_width
    for col in range(block_width):
        resnum = start_resnum + col
        if resnum % 10 != 0:
            continue
        ticks[col] = "|"
        digits = str(resnum)
        start = col - len(digits) + 1
        for i, d in enumerate(digits):
            pos = start + i
            if 0 <= pos < block_width:
                numbers[pos] = d
    return "".join(numbers), "".join(ticks)


def _one_letter(resname: str) -> str:
    return AAMAP.get(resname, "X")


def _find_structure(structures: list[LabeledStructure], label: str) -> LabeledStructure:
    for s in structures:
        if s.label == label:
            return s
    raise ValueError(f"structure not found: {label!r} (labels resolved via --indir: {[s.label for s in structures]})")


def _find_chain(structure: LabeledStructure, chain_id: str) -> ChainSequence | None:
    for c in structure.chains:
        if c.chain_id == chain_id:
            return c
    return None


def _format_range(start: int, end: int) -> str:
    return f"{start}-{end}" if start != end else str(start)


def _format_resnum_ranges(resnums: set[int] | list[int]) -> str:
    """残基番号の集合を'12-18, 25, 40-42'のような連続範囲表記に整形する。"""
    ordered = sorted(resnums)
    if not ordered:
        return ""
    ranges = []
    start = prev = ordered[0]
    for n in ordered[1:]:
        if n == prev + 1:
            prev = n
            continue
        ranges.append(_format_range(start, prev))
        start = prev = n
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ", ".join(ranges)


def _format_gap_note(ref_chain: ChainSequence | None, other_chain: ChainSequence | None) -> str | None:
    """基準チェーンと対象チェーンの残基番号集合を比較し、欠損(片方にしかない)範囲を整形する。

    PDB構造は電子密度が見えない領域(ループ等)が欠落することが多いため、置換とは別に
    可視化する。どちらのチェーンも見つからない場合は`None`を返す。
    """
    if ref_chain is None or other_chain is None:
        return None
    ref_set, other_set = set(ref_chain.resnums), set(other_chain.resnums)
    missing_in_other = ref_set - other_set
    extra_in_other = other_set - ref_set
    notes = []
    if missing_in_other:
        notes.append(f"reference only (missing in target): {_format_resnum_ranges(missing_in_other)}")
    if extra_in_other:
        notes.append(f"target only (absent in reference): {_format_resnum_ranges(extra_in_other)}")
    if not notes:
        return None
    return "; ".join(notes)


def format_mutation_report(structures: list[LabeledStructure], reference: str) -> str:
    """`reference`を基準に、他の構造との残基置換一覧を整形する。

    `reference`は次のいずれか:
    - \"ラベル:チェーンID\"(例: `P24941_AF:A`): `--indir`等で読み込んだ構造の1チェーンを
      基準にする。残基番号ベースの対応付けのみを用いる(構造間でPDBの残基番号が
      揃っている前提。`pf align-view --method number`/`structfit`と同じ前提、
      `structcompare`参照)。
    - アミノ酸配列(1文字表記、コロンを含まない): 構造を伴わない任意配列
      (UniProt正規配列等)を基準にする。この場合は`seqalign`によるグローバル
      配列アラインメントを用いるため、残基番号が揃っていない構造間でも比較できる。
    """
    if ":" in reference:
        return _format_mutation_report_vs_structure(structures, reference)
    return _format_mutation_report_vs_sequence(structures, reference)


def _format_mutation_report_vs_structure(structures: list[LabeledStructure], reference: str) -> str:
    ref_label, _, ref_chain_id = reference.partition(":")
    ref_structure = _find_structure(structures, ref_label)
    ref_atoms = ref_structure.atoms.select(f"protein and chain {ref_chain_id}")
    if ref_atoms is None:
        raise ValueError(f"chain not found: {reference!r}")
    ref_chain = _find_chain(ref_structure, ref_chain_id)

    lines = [f"reference: {reference}"]
    for s in structures:
        if s.label == ref_label:
            continue
        report = find_substitutions(ref_atoms, s.atoms)
        if not report.matched:
            lines.append(f"  {s.label}: no residue-number-based correspondence found (numbering scheme may differ)")
            continue
        if not report.substitutions:
            lines.append(f"  {s.label}: no substitutions (seqid={report.seqid:.1f}%, overlap={report.overlap:.1f}%)")
        else:
            subs_str = ", ".join(
                f"{_one_letter(sub.resname_a)}{sub.resnum}{_one_letter(sub.resname_b)}"
                for sub in report.substitutions
            )
            lines.append(
                f"  {s.label}: {len(report.substitutions)} substitution(s)"
                f" (seqid={report.seqid:.1f}%, overlap={report.overlap:.1f}%): {subs_str}"
            )
        gap_note = _format_gap_note(ref_chain, _find_chain(s, report.chain_id_b))
        if gap_note:
            lines.append(f"    gaps: {gap_note}")
    return "\n".join(lines) + "\n"


def _format_mutation_report_vs_sequence(structures: list[LabeledStructure], reference: str) -> str:
    ref_sequence = reference.upper()
    if not _SEQUENCE_PATTERN.match(ref_sequence):
        raise ValueError(
            "--reference must be either 'label:chain_id' (e.g. P24941_AF:A) "
            f"or an amino acid sequence (one-letter code): {reference!r}"
        )

    lines = [f"reference sequence: {len(ref_sequence)} residues"]
    for s in structures:
        if not s.chains:
            lines.append(f"  {s.label}: no protein chains found")
            continue
        for c in s.chains:
            result = align_to_reference(ref_sequence, c.sequence, c.resnums)
            label = f"{s.label}:{c.chain_id}"
            stats = f"identity={result.identity:.1f}%, coverage={result.coverage:.1f}%"
            if not result.substitutions:
                lines.append(f"  {label}: no substitutions ({stats})")
            else:
                # 基準配列内の位置(ref_pos)と構造側の実際の残基番号(query_resnum)は
                # 一致するとは限らない(基準配列がUniProt正規配列と一致しない場合等)ため、
                # 両方を明示する。
                subs_str = ", ".join(
                    f"{sub.ref_aa}{sub.ref_pos}{sub.query_aa}(structure resnum={sub.query_resnum})"
                    for sub in result.substitutions
                )
                lines.append(f"  {label}: {len(result.substitutions)} substitution(s) ({stats}): {subs_str}")
            deletions = [g for g in result.gaps if g.kind == "deletion"]
            if deletions:
                del_str = ", ".join(
                    f"ref {_format_range(g.ref_start, g.ref_end)} ({g.length} residue(s))" for g in deletions
                )
                lines.append(f"    gaps: {del_str}")
            insertions = [g for g in result.gaps if g.kind == "insertion"]
            if insertions:
                ins_str = ", ".join(f"{g.length} residue(s)" for g in insertions)
                lines.append(f"    not in reference sequence: {ins_str}")
    return "\n".join(lines) + "\n"


def build_report(structures: list[LabeledStructure], reference: str | None) -> str:
    """FASTA・pairwise identity・整列表示・(reference指定時)置換一覧をまとめたレポートを組み立てる。"""
    parts = [
        "== Sequences (FASTA, observed residues only) ==",
        format_fasta(structures),
        "== Pairwise identity ==",
        format_identity_matrix(structures),
        "== Alignment (by residue number) ==",
        format_alignment_block(structures),
    ]
    if reference:
        parts += ["== Substitutions relative to reference ==", format_mutation_report(structures, reference)]
    return "\n".join(parts)

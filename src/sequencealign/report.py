"""複数構造間の蛋白配列比較レポート(FASTA・pairwise identity・基準配列に対する置換一覧)。"""

import re
import string
from dataclasses import dataclass
from pathlib import Path

from prody.atomic.atomic import AAMAP, Atomic

from core.logging_utils import get_logger
from seqalign import AlignmentResult, align_to_reference
from seqextract import ChainSequence, get_chain_sequences
from structcompare import find_substitutions
from structio import parse_fasta, parse_structure

logger = get_logger(__name__)

_FASTA_WIDTH = 60
DEFAULT_ALIGN_WIDTH = 100
# 標準20種 + 曖昧/非標準コード(Asx/Glx/Xle/Sec/Pyl/不明)。matchChains側のAAMAPと合わせる。
_SEQUENCE_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZJUO]+$")


@dataclass
class LabeledStructure:
    """入力トークンから解決したファイル名(拡張子抜き)をラベルとする構造。

    `.fasta`から読み込んだ場合、`atoms`は`None`になる(3次元構造を持たないため)。
    pairwise identity(`format_identity_matrix`)や構造基準の置換検出
    (`format_mutation_report`の`label:chain_id`形式)など`atoms`を要する処理は、
    `atoms`が`None`の構造を対象から除外する。
    """

    label: str
    atoms: Atomic | None
    chains: list[ChainSequence]


def load_labeled_structures(paths: list[Path]) -> list[LabeledStructure]:
    """構造ファイル(PDB/CIF)またはFASTAファイルを読み込み、ファイル名(拡張子抜き)を
    ラベルとして付与する。
    """
    structures = []
    for path in paths:
        if path.suffix.lower() == ".fasta":
            logger.info("Loading sequence(s) from %s ...", path)
            chains = _chains_from_fasta(path)
            logger.info("  -> %d sequence(s): %s", len(chains), [c.chain_id for c in chains])
            structures.append(LabeledStructure(label=path.stem, atoms=None, chains=chains))
        else:
            logger.info("Loading structure from %s ...", path)
            atoms = parse_structure(path)
            chains = get_chain_sequences(atoms)
            logger.info("  -> %d protein chain(s): %s", len(chains), [c.chain_id for c in chains])
            structures.append(LabeledStructure(label=path.stem, atoms=atoms, chains=chains))
    return structures


def _chains_from_fasta(path: Path) -> list[ChainSequence]:
    """FASTAの各レコードを、A/B/C...と順にチェーンIDを振った`ChainSequence`に変換する。

    FASTAには残基番号情報がないため、1残基目をresnum=1として連番を振る
    (`format_alignment_block`のsequenceタイプreferenceと同じ前提)。
    """
    records = parse_fasta(path)
    chains = []
    for i, (_, sequence) in enumerate(records):
        chain_id = string.ascii_uppercase[i] if i < len(string.ascii_uppercase) else str(i + 1)
        chains.append(ChainSequence(chain_id=chain_id, sequence=sequence, resnums=list(range(1, len(sequence) + 1))))
    return chains


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


def _pairwise_grid_entries(structures: list[LabeledStructure]) -> list[tuple[str, ChainSequence]]:
    return [(f"{s.label}:{c.chain_id}", c) for s in structures for c in s.chains]


def _pairwise_alignment_grid(
    structures: list[LabeledStructure],
) -> tuple[list[str], list[list[AlignmentResult | None]]]:
    """全チェーンの組み合わせ(順序あり、対角除く)についてペアワイズグローバルアラインメント
    (`seqalign.align_to_reference`)を計算する。`grid[i][j]`は行iを基準配列(ref)、列jを
    比較対象(query)とした`AlignmentResult`(対角は`None`)。identityはほぼ対称だが、
    coverageは基準側の長さに対する割合のため方向で値が異なるため、両方向を計算する。
    """
    entries = _pairwise_grid_entries(structures)
    labels = [label for label, _ in entries]
    n = len(labels)
    grid: list[list[AlignmentResult | None]] = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            ref_chain, query_chain = entries[i][1], entries[j][1]
            grid[i][j] = align_to_reference(ref_chain.sequence, query_chain.sequence, query_chain.resnums)
    return labels, grid


def _render_grid_table(labels: list[str], cells: list[list[str]]) -> str:
    """`labels`を行・列ラベルとするN×Nのグリッド表を整形する(`cells[i][j]`は行i・列jの値)。"""
    row_label_width = max(len(label) for label in labels)
    col_width = max([row_label_width] + [len(v) for row in cells for v in row])

    header = " " * (row_label_width + 2) + "".join(label.rjust(col_width + 1) for label in labels)
    lines = [header]
    for i, label in enumerate(labels):
        row_cells = "".join(cells[i][j].rjust(col_width + 1) for j in range(len(labels)))
        lines.append(f"{label.ljust(row_label_width)}  {row_cells}")
    return "\n".join(lines) + "\n"


def format_identity_matrix(structures: list[LabeledStructure]) -> str:
    """全チェーンの組み合わせについて、%identityをN×Nのグリッド表として一覧化する
    (対角は自身との比較のため`-`)。

    構造(`atoms`)の有無によらず、全チェーンの配列同士をペアワイズグローバルアラインメント
    (`seqalign.align_to_reference`)し、アラインメントされた位置における%identityをセルの
    値とする。同一構造内の複数チェーン同士の組(例: ホモ二量体のチェーンA/C)も、FASTA由来の
    チェーン(atoms無し)を含む組も対象にする。identityは基準・対象どちらを基準にしても
    ほぼ同じ値になるため、対角より上のセルのみ計算し下は複製する(`format_coverage_matrix`とは
    異なり非対称ではない)。
    """
    entries = _pairwise_grid_entries(structures)
    if len(entries) < 2:
        return "(fewer than two protein chains found)\n"

    labels = [label for label, _ in entries]
    n = len(labels)
    cells = [["-"] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            result = align_to_reference(entries[i][1].sequence, entries[j][1].sequence, entries[j][1].resnums)
            value = f"{result.identity:.1f}"
            cells[i][j] = value
            cells[j][i] = value
    return _render_grid_table(labels, cells)


def format_coverage_matrix(structures: list[LabeledStructure]) -> str:
    """全チェーンの組み合わせについて、%coverage(アラインメントされた割合)をN×Nの
    グリッド表として一覧化する(対角は自身との比較のため`-`)。

    セル`(行, 列)`は「行のチェーンを基準配列としたとき、アラインメントで列のチェーンと
    対応した割合」を表す(基準配列の長さが分母のため、行と列を入れ替えると値が変わる
    非対称な表。例: 短いドメインのみの構造を基準にすると高coverage、全長配列を基準にすると
    そのドメインの分しかカバーされないため低coverageになる)。
    """
    labels, grid = _pairwise_alignment_grid(structures)
    if len(labels) < 2:
        return "(fewer than two protein chains found)\n"

    n = len(labels)
    cells = [[grid[i][j].coverage if grid[i][j] else None for j in range(n)] for i in range(n)]
    formatted = [[f"{v:.1f}" if v is not None else "-" for v in row] for row in cells]
    return _render_grid_table(labels, formatted)


def format_identity_coverage_matrix(structures: list[LabeledStructure]) -> str:
    """全チェーンの組み合わせについて、`identity/coverage`をN×Nのグリッド表として
    一覧化する(対角は自身との比較のため`-`)。`format_identity_matrix`と
    `format_coverage_matrix`を1つの表に統合した形式(coverageの非対称性については
    `format_coverage_matrix`のdocstring参照)。
    """
    labels, grid = _pairwise_alignment_grid(structures)
    if len(labels) < 2:
        return "(fewer than two protein chains found)\n"

    n = len(labels)
    cells = [
        [f"{grid[i][j].identity:.1f}/{grid[i][j].coverage:.1f}" if grid[i][j] else "-" for j in range(n)]
        for i in range(n)
    ]
    return _render_grid_table(labels, cells)


def format_alignment_block(structures: list[LabeledStructure], width: int = DEFAULT_ALIGN_WIDTH) -> str:
    """全構造・全蛋白チェーンの配列を、残基番号を共通の軸として横並びに整列表示する
    (`width`残基ごとに折り返す)。

    配列アラインメントは行わず、残基番号が一致する列に同じアミノ酸が並ぶ前提で
    並べる(構造間でPDBの残基番号が揃っている前提。`pf align-view --method number`
    と同じ前提)。観測されていない残基は`-`で埋める。番号体系が揃っていない構造
    (例: 全長のUniProt配列と、その一部ドメインのみを含む結晶構造)を混在させると
    無意味な結果になる点に注意(その場合は`format_alignment_block_by_sequence`を使う)。

    基準配列(UniProt正規配列等)を加えたい場合は、その配列を`.fasta`ファイルとして
    入力トークンに含めればよい(`load_labeled_structures`が1残基目をresnum=1として
    読み込むため、他の行と同じ軸で並ぶ。特別な「reference」行の仕組みは持たない)。
    """
    entries: list[tuple[str, str, list[int]]] = [
        (f"{s.label}:{c.chain_id}", c.sequence, c.resnums) for s in structures for c in s.chains
    ]
    if not entries:
        return "(no protein chains found)\n"

    all_resnums = {r for _, _, resnums in entries for r in resnums}
    min_resnum, max_resnum = min(all_resnums), max(all_resnums)

    rows = []
    for label, sequence, resnums in entries:
        seq_by_resnum = dict(zip(resnums, sequence))
        padded = "".join(seq_by_resnum.get(r, "-") for r in range(min_resnum, max_resnum + 1))
        rows.append((label, padded))

    return _render_alignment_blocks(rows, width, start_pos=min_resnum)


def format_alignment_block_by_sequence(structures: list[LabeledStructure], width: int = DEFAULT_ALIGN_WIDTH) -> str:
    """全構造・全蛋白チェーンの配列を、配列の相同性(ペアワイズグローバルアラインメント)に
    基づいて位置を揃えて横並びに整列表示する(`width`残基ごとに折り返す)。

    `format_alignment_block`(残基番号ベース)と異なり、構造間でPDBの残基番号が揃っていない
    場合でも正しい位置に並べられる(例: 全長のUniProt配列と、その一部ドメインのみを含む
    結晶構造の組み合わせ)。基準は先頭の構造の先頭チェーン(`pf align-view --method align`が
    先頭構造を基準にするのと同じ考え方)とし、他の全チェーンをこの基準配列に対して個別に
    ペアワイズアラインメントする(`seqalign.align_to_reference`、Biopython
    `PairwiseAligner`によるグローバルアラインメント)ため、複数配列同時アラインメント(MSA)
    ではない点に注意。異なる蛋白の配列を混在させると無意味な結果になる点は
    `format_alignment_block`と同じ。
    """
    all_chains = [(f"{s.label}:{c.chain_id}", c) for s in structures for c in s.chains]
    if not all_chains:
        return "(no protein chains found)\n"

    ref_label, ref_chain = all_chains[0]
    ref_sequence = ref_chain.sequence
    ref_length = len(ref_sequence)

    rows = [(ref_label, ref_sequence)]
    for label, chain in all_chains[1:]:
        result = align_to_reference(ref_sequence, chain.sequence, chain.resnums)
        padded = "".join(result.query_by_ref_pos.get(pos, "-") for pos in range(1, ref_length + 1))
        rows.append((label, padded))

    return _render_alignment_blocks(rows, width, start_pos=1)


def _render_alignment_blocks(rows: list[tuple[str, str]], width: int, start_pos: int) -> str:
    """`(label, padded_sequence)`の行リストを、`width`残基ごとにルーラー付きで折り返す。

    全行の`padded_sequence`は同じ長さで、列`i`(0始まり)が軸上の位置`start_pos + i`に
    対応する前提(`format_alignment_block`/`format_alignment_block_by_sequence`共通)。
    """
    label_width = max(len(label) for label, _ in rows)
    total_length = len(rows[0][1])
    blocks = []
    for block_start in range(0, total_length, width):
        block_end = min(block_start + width, total_length)
        block_first_pos = start_pos + block_start
        number_line, tick_line = _format_ruler(block_first_pos, block_end - block_start)
        indent = " " * (label_width + 2)
        block_lines = [indent + number_line, indent + tick_line]
        block_lines += [f"{label.ljust(label_width)}  {seq[block_start:block_end]}" for label, seq in rows]
        blocks.append("\n".join(block_lines))
    return "\n\n".join(blocks) + "\n"


def _format_ruler(start_resnum: int, block_width: int) -> tuple[str, str]:
    """10残基ごとに残基番号とその位置を示す目盛り(2行: 数字の行、`|`の行)を作る。

    数字はその残基番号の列で右端が揃うように配置する(例: 残基120の場合、'0'が
    resnum=120の列に来る)。ただしブロック左端に近く数字全体が右揃えでは収まらない
    目盛り(通常は各ブロック最初の目盛りのみ)は、上位の桁が欠けて末尾の'0'だけが
    見える(実際の値を誤読させる)ことを避けるため、ブロック左端(列0)に寄せて
    数字全体を表示する(`|`自体は本来の列のまま動かさない)。隣接する目盛り同士は
    10列以上離れており数字は最大4桁のため、この寄せによる重なりは生じない。
    """
    numbers = [" "] * block_width
    ticks = [" "] * block_width
    for col in range(block_width):
        resnum = start_resnum + col
        if resnum % 10 != 0:
            continue
        ticks[col] = "|"
        digits = str(resnum)
        start = max(0, col - len(digits) + 1)
        for i, d in enumerate(digits):
            pos = start + i
            if pos < block_width:
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
    if ref_structure.atoms is None:
        raise ValueError(
            f"reference {reference!r} has no atomic structure (loaded from FASTA); "
            "use an amino acid sequence as --reference instead"
        )
    ref_atoms = ref_structure.atoms.select(f"protein and chain {ref_chain_id}")
    if ref_atoms is None:
        raise ValueError(f"chain not found: {reference!r}")
    ref_chain = _find_chain(ref_structure, ref_chain_id)

    lines = [f"reference: {reference}", _wrap(ref_chain.sequence) if ref_chain else ""]
    for s in structures:
        if s.label == ref_label:
            continue
        if s.atoms is None:
            lines.append(f"  {s.label}: no atomic structure (loaded from FASTA); residue-number-based comparison not available")
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

    lines = [f"reference sequence: {len(ref_sequence)} residues", _wrap(ref_sequence)]
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


def build_report(
    structures: list[LabeledStructure],
    align_width: int = DEFAULT_ALIGN_WIDTH,
    method: str = "align",
    identity_format: str = "combined",
) -> str:
    """Pairwise identity・整列表示をまとめたレポートを組み立てる。

    `method`(`align`/`number`)は整列表示の方式を選ぶ。`align`(既定)は配列の相同性に基づく
    ペアワイズアラインメントベース(`format_alignment_block_by_sequence`)で、番号体系が
    揃っていない構造の組み合わせでも正しく並ぶ。`number`は残基番号ベース
    (`format_alignment_block`)で、構造間でPDBの残基番号が既に揃っている(同じUniProt番号体系
    等)ことが分かっている場合にのみ使う。

    `identity_format`(`combined`/`separate`)はPairwise identityセクションの表示形式を選ぶ。
    `combined`(既定)はidentity/coverageを1つの表(セルは`identity/coverage`)にまとめる。
    `separate`はidentity表・coverage表を別々のセクションとして出力する
    (どちらの形式にするか、実際の出力を見て決めたいというユーザーの要望により両方実装した)。
    """
    if method == "number":
        alignment_heading = "== Alignment (by residue number) =="
        alignment = format_alignment_block(structures, width=align_width)
    else:
        alignment_heading = "== Alignment (sequence-aligned) =="
        alignment = format_alignment_block_by_sequence(structures, width=align_width)

    if identity_format == "separate":
        parts = [
            "== Pairwise identity ==",
            format_identity_matrix(structures),
            "== Coverage ==",
            format_coverage_matrix(structures),
        ]
    else:
        parts = [
            "== Pairwise identity/coverage ==",
            format_identity_coverage_matrix(structures),
        ]
    parts += [alignment_heading, alignment]
    return "\n".join(parts)

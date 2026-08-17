"""検出済みのポケット(list[Pocket])を、周辺残基を強調表示した状態でPyMOLで開く。"""

from pathlib import Path

from pocket import Pocket, PocketResidue
from pymolrun import run_pymol_script

from core.logging_utils import get_logger

logger = get_logger(__name__)

# ポケットを見分けやすいよう、スコア順に割り当てる配色。
_POCKET_COLORS = ["red", "orange", "yellow", "green", "cyan", "blue", "magenta", "purple"]


def _residue_selection(residues: list[PocketResidue]) -> str:
    """残基一覧を、チェーンごとにまとめたPyMOL選択式(`(chain A and resi 1+2+..) or (chain B and resi ..)`)にする。"""
    by_chain: dict[str, list[int]] = {}
    for residue in residues:
        by_chain.setdefault(residue.chain_id, []).append(residue.resnum)
    return " or ".join(
        f"(chain {chain} and resi {'+'.join(str(resnum) for resnum in sorted(resnums))})"
        for chain, resnums in by_chain.items()
    )


def build_pocket_view_script(structure_path: Path, pockets: list[Pocket], top_n: int | None = None) -> str:
    """構造を読み込み、ポケットごとに周辺残基を配色・強調表示するPyMOLスクリプト(.pml)本文を組み立てる。

    `pockets`は`report.read_pockets_json()`のスコア降順の出力を想定。`top_n`指定時は上位N件のみ強調する。
    """
    if top_n is not None:
        pockets = pockets[:top_n]

    object_name = Path(structure_path).stem
    lines = [f"load {Path(structure_path).resolve()}, {object_name}"]
    lines += ["hide everything", f"show cartoon, {object_name}", f"color gray80, {object_name}"]

    colors = (_POCKET_COLORS * (len(pockets) // len(_POCKET_COLORS) + 1))[: len(pockets)]
    highlighted_names = []
    for pocket, color in zip(pockets, colors):
        if not pocket.residues:
            logger.warning("Pocket %d has no residues; skipping in view", pocket.pocket_id)
            continue
        name = f"pocket_{pocket.pocket_id}"
        highlighted_names.append(name)
        lines += [
            f"select {name}, {object_name} and ({_residue_selection(pocket.residues)})",
            f"color {color}, {name}",
            f"show sticks, {name} and sidechain",
            f"show surface, {name}",
            f"set transparency, 0.3, {name}",
            f'print("view-pocket: pocket {pocket.pocket_id} score={pocket.score:.3f} '
            f'druggability={pocket.druggability_score:.3f} ({len(pocket.residues)} residues)")',
        ]

    lines.append(f"zoom ({' or '.join(highlighted_names)})" if highlighted_names else f"zoom {object_name}")
    lines.append("deselect")
    return "\n".join(lines) + "\n"


def launch_pocket_view(
    structure_path: Path, pockets: list[Pocket], pymol_env: str = "pymol", top_n: int | None = None
) -> None:
    """構造とポケット一覧をPyMOLで開き、ポケットごとに周辺残基を強調表示する。"""
    logger.info("Building PyMOL script for %d pocket(s) ...", len(pockets) if top_n is None else top_n)
    script = build_pocket_view_script(structure_path, pockets, top_n=top_n)
    run_pymol_script(script, pymol_env=pymol_env)

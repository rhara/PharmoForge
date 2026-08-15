"""構造ファイルから指定チェーンの抽出・水分子の除去を行う。"""

from pathlib import Path

from core.logging_utils import get_logger
from structio import parse_structure, write_structure

logger = get_logger(__name__)


def extract_structure(
    input_path: Path,
    output_path: Path,
    chains: list[str] | None = None,
    remove_water: bool = False,
) -> None:
    """`input_path`の構造から指定チェーンを抽出し(`chains`省略時は全チェーン)、
    `remove_water=True`なら水分子を除いて`output_path`に保存する。

    入出力とも拡張子(`.pdb`/`.cif`)で形式を自動判別するため、入力と出力で異なる
    形式を指定してもよい。
    """
    input_path = Path(input_path)
    logger.info("Loading structure from %s ...", input_path)
    structure = parse_structure(input_path)

    selection_parts = []
    if chains:
        selection_parts.append("chain " + " ".join(chains))
    if remove_water:
        selection_parts.append("not water")
    selection = " and ".join(selection_parts) if selection_parts else "all"

    logger.info("Selecting atoms (%s) ...", selection)
    selected = structure.select(selection)
    if selected is None:
        raise ValueError(f"selection matched no atoms: {selection!r}")

    logger.info(
        "Selected %d / %d atoms (chains=%s)",
        len(selected), structure.numAtoms(), sorted(set(selected.getChids())),
    )
    write_structure(selected, Path(output_path))
    logger.info("Done: saved extracted structure to %s", output_path)

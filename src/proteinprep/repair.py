"""PDBFixerによる蛋白構造の基本修復(欠損原子補完・プロトン化)。"""

from pathlib import Path

from openmm.app import PDBFile
from pdbfixer import PDBFixer

from core.logging_utils import get_logger

logger = get_logger(__name__)


def repair_structure(input_path: Path, output_path: Path, ph: float | None = None) -> None:
    """蛋白構造の欠損残基・欠損原子を補完し、output_pathにPDB形式で保存する。

    ph を指定した場合(MD用途)、その pH で水素原子を付加する。
    ph=None(既定、docking用途)の場合は水素原子を付加しない。
    リガンドや水などのヘテロ原子は変更・除去しない(分類・整形は別機能で扱う)。
    入力はPDBFixerが受け付ける形式(PDB/CIF)。
    """
    logger.info("Loading structure from %s ...", input_path)
    fixer = PDBFixer(filename=str(input_path))

    logger.info("Finding missing residues and atoms ...")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    logger.info("Adding missing atoms ...")
    fixer.addMissingAtoms()

    if ph is not None:
        logger.info("Adding hydrogens at pH %.2f ...", ph)
        fixer.addMissingHydrogens(ph)
    else:
        logger.info("Skipping hydrogen addition (docking mode)")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    logger.info("Done: saved repaired structure to %s", output_path)

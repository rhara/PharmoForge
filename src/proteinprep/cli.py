from pathlib import Path

import click

from .repair import repair_structure


@click.command("prep-protein")
@click.argument("input_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    "output_path",
    required=True,
    type=click.Path(path_type=Path),
    help="出力PDBファイルパス",
)
@click.option(
    "--mode",
    type=click.Choice(["dock", "md"]),
    default="dock",
    show_default=True,
    help="dock=水素原子を付加しない、md=指定pHでプロトン化する",
)
@click.option(
    "--ph",
    type=float,
    default=7.0,
    show_default=True,
    help="--mode md 時のプロトン化pH",
)
def prep_protein_cmd(input_path: Path, output_path: Path, mode: str, ph: float):
    """蛋白構造(PDB/CIF)の欠損原子を補完し、必要に応じてプロトン化する(PDBFixer)。

    \b
    例:
      pf prep-protein data/9csk.cif --output data/9csk_dock.pdb --mode dock
      pf prep-protein data/P61626.cif --output data/P61626_md.pdb --mode md --ph 7.4
    """
    ph_value = ph if mode == "md" else None
    repair_structure(input_path, output_path, ph=ph_value)

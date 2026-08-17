from pathlib import Path

import click

from pocket import run_fpocket

from .report import format_pockets_table, read_pockets_json, write_pockets_json, write_pockets_table
from .view import launch_pocket_view


@click.command("find-pocket")
@click.argument("structure_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--outdir",
    "-o",
    "outdir",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for fpocket's raw output and the pocket list (pockets.json)",
)
@click.option(
    "--top",
    "top_n",
    default=None,
    type=int,
    help="Only write the top N pockets by score to pockets.json (all pockets if omitted)",
)
def find_pocket_cmd(structure_path: Path, outdir: Path, top_n: int | None):
    """Detect candidate ligand-binding pockets in a protein structure (PDB/CIF) and list the residues lining each pocket.

    OUTDIR keeps the pocket list (`pockets.json`) alongside fpocket's raw output
    (`<structure name>_out/`, including PyMOL/VMD visualization scripts).

    \b
    Examples:
      pf find-pocket data/cdk2/P24941_AF.cif --outdir data/cdk2/P24941_AF_pockets
      pf find-pocket data/cdk2/P24941_AF.cif --outdir data/cdk2/P24941_AF_pockets --top 3
    """
    pockets = run_fpocket(structure_path, outdir)
    if top_n is not None:
        pockets = pockets[:top_n]
    write_pockets_json(structure_path.name, pockets, outdir / "pockets.json")


@click.command("view-pocket")
@click.argument("structure_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--pockets",
    "-p",
    "pockets_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to the pockets.json produced by `pf find-pocket`",
)
@click.option(
    "--top",
    "top_n",
    default=None,
    type=int,
    help="Highlight only the top N pockets by score (all pockets in the file if omitted)",
)
@click.option(
    "--pymol-env",
    default="pymol",
    show_default=True,
    help="Name of the dedicated conda/mamba environment PyMOL is installed in",
)
def view_pocket_cmd(structure_path: Path, pockets_path: Path, top_n: int | None, pymol_env: str):
    """Open a structure (PDB/CIF) in PyMOL with the residues of each detected pocket highlighted.

    STRUCTURE_PATH should be the same structure `pf find-pocket` was run on (or an equivalent
    one with matching chain IDs/residue numbers); --pockets points to its pockets.json output.

    \b
    Examples:
      pf view-pocket data/cdk2/P24941_AF.cif --pockets data/cdk2/P24941_AF_pockets/pockets.json
      pf view-pocket data/cdk2/P24941_AF.cif --pockets data/cdk2/P24941_AF_pockets/pockets.json --top 3
    """
    pockets = read_pockets_json(pockets_path)
    launch_pocket_view(structure_path, pockets, pymol_env=pymol_env, top_n=top_n)


@click.command("show-pocket")
@click.argument("pockets_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    "output_path",
    default=None,
    type=click.Path(path_type=Path),
    help="Output TSV file path (prints to stdout if omitted)",
)
def show_pocket_cmd(pockets_path: Path, output_path: Path | None):
    """Print the pockets.json produced by `pf find-pocket` as a table (TSV, one row per pocket).

    The `residues` column packs each pocket's residues as a comma-separated `<chain>:<resnum>` list.

    \b
    Examples:
      pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json
      pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json --output data/cdk2/P24941_AF_pockets/pockets.tsv
      pf show-pocket data/cdk2/P24941_AF_pockets/pockets.json | column -t
    """
    pockets = read_pockets_json(pockets_path)
    if output_path is None:
        click.echo(format_pockets_table(pockets), nl=False)
    else:
        write_pockets_table(pockets, output_path)

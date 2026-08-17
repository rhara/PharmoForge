from pathlib import Path

import click

from .extract import extract_structure


@click.command("protein-extract")
@click.argument("input_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--chains",
    default=None,
    help="Comma-separated chain IDs to extract (auth_asym_id, the ID shown in PyMOL/RCSB's website; e.g. A,B). All chains if omitted.",
)
@click.option(
    "--remove-water",
    is_flag=True,
    default=False,
    help="Remove water molecules (HOH etc.)",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    required=True,
    type=click.Path(path_type=Path),
    help="Output structure file path (format determined by the .pdb/.cif extension)",
)
def protein_extract_cmd(input_path: Path, chains: str | None, remove_water: bool, output_path: Path):
    """Extract the given chains from a structure file (PDB/CIF), optionally removing water molecules.

    Both input and output formats (`.pdb`/`.cif`) are auto-detected from the extension, so this also
    converts between formats.

    \b
    Examples:
      pf protein-extract data/cdk2/2CCH.cif --chains=A --remove-water --output data/cdk2/2CCH_main.cif
      pf protein-extract data/tyk2/6NZP.cif --remove-water --output data/tyk2/6NZP_nowater.pdb
    """
    chain_list = [c.strip() for c in chains.split(",") if c.strip()] if chains else None
    extract_structure(input_path, output_path, chains=chain_list, remove_water=remove_water)

from pathlib import Path

import click

from idmap.identifiers import entry_name_to_accession, looks_like_uniprot_accession
from uniprot import fetch_protein_info

from .report import format_protein_info_json, write_protein_info_json


@click.command("protein-info")
@click.argument("identifier")
@click.option(
    "--output",
    "-o",
    "output_path",
    default=None,
    type=click.Path(path_type=Path),
    help="Output JSON file path (prints JSON to stdout if omitted)",
)
def protein_info_cmd(identifier: str, output_path: Path | None):
    """Fetch drug-discovery-relevant protein info (structural biology/medicinal chemistry) from UniProt and save it as JSON.

    IDENTIFIER is a UniProt entry name (e.g. EGFR_HUMAN) or accession (e.g. P00533).

    \b
    Examples:
      pf protein-info EGFR_HUMAN --output data/egfr_info.json
      pf protein-info P00533 --output data/egfr_info.json
      pf protein-info EGFR_HUMAN
      pf protein-info EGFR_HUMAN | jq .pdb_structures
    """
    identifier = identifier.strip()
    accession = (
        identifier if looks_like_uniprot_accession(identifier) else entry_name_to_accession(identifier)
    )
    info = fetch_protein_info(accession)
    if output_path is None:
        click.echo(format_protein_info_json(info))
    else:
        write_protein_info_json(info, output_path)

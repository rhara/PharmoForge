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
    help="出力JSONファイルパス(省略時は標準出力にJSONを出力する)",
)
def protein_info_cmd(identifier: str, output_path: Path | None):
    """UniProtから創薬(構造生物学・メディシナルケミストリー)向けの蛋白情報を取得しJSONで保存する。

    IDENTIFIERにはUniProt entry name(例: EGFR_HUMAN)またはaccession(例: P00533)を指定する。

    \b
    例:
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

from pathlib import Path

import click

from idmap.identifiers import entry_name_to_accession, looks_like_uniprot_accession
from uniprot import fetch_protein_info

from .report import write_protein_info_json


@click.command("protein-info")
@click.argument("identifier")
@click.option(
    "--output",
    "-o",
    "output_path",
    required=True,
    type=click.Path(path_type=Path),
    help="出力JSONファイルパス",
)
def protein_info_cmd(identifier: str, output_path: Path):
    """UniProtから創薬(構造生物学・メディシナルケミストリー)向けの蛋白情報を取得しJSONで保存する。

    IDENTIFIERにはUniProt entry name(例: EGFR_HUMAN)またはaccession(例: P00533)を指定する。

    \b
    例:
      pf protein-info EGFR_HUMAN --output data/egfr_info.json
      pf protein-info P00533 --output data/egfr_info.json
    """
    identifier = identifier.strip()
    accession = (
        identifier if looks_like_uniprot_accession(identifier) else entry_name_to_accession(identifier)
    )
    info = fetch_protein_info(accession)
    write_protein_info_json(info, output_path)

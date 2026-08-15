from pathlib import Path

import click

from afdb import fetch_structure as fetch_af_structure
from afdb import fetch_structures as fetch_af_structures
from chembl import fetch_activities
from idmap import resolve_chembl_target_id, resolve_uniprot_accession
from rcsb import fetch_structure, fetch_structures

from . import activities

SUPPORTED_TYPES = ("activities", "structure", "structures")


@click.command("fetch")
@click.argument("spec")
@click.option(
    "--af",
    is_flag=True,
    default=False,
    help="structure=/structures=の取得元をAlphaFold DBにする(既定はRCSB PDB)。"
    "識別子はUniProt entry name/accessionを指定する。",
)
@click.option(
    "--output",
    "-o",
    "output",
    required=True,
    type=click.Path(path_type=Path),
    help="出力ファイルパス(structures=ではディレクトリ)",
)
@click.option(
    "--type",
    "fmt",
    type=click.Choice(["cif", "pdb"], case_sensitive=False),
    default=None,
    help="構造ファイルのフォーマット(structure/structures用)。structure=単体では省略時にoutputの拡張子から推定する。",
)
def fetch_cmd(spec: str, af: bool, output: Path, fmt: str | None):
    """データを取得する。SPECは <データ種別>=<識別子> の形式。

    \b
    例:
      pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
      pf fetch structure=9CSK --type=cif --output data/9csk.cif
      pf fetch structures=6P8F,7SJ3,9CSK --type pdb --output data
      pf fetch structure=P61626 --af --type=cif --output data/P61626.cif
      pf fetch structure=TYK2_HUMAN --af --type=cif --output data/tyk2/TYK2_HUMAN_af.cif
      pf fetch structures=P61626,CDK4_HUMAN --af --type pdb --output data
    """
    if "=" not in spec:
        raise click.UsageError(
            f"SPECは <データ種別>=<識別子> の形式で指定してください(例: activities=CDK4_HUMAN): {spec!r}"
        )
    data_type, _, value = spec.partition("=")
    data_type = data_type.strip().lower()
    value = value.strip()

    if af and data_type not in ("structure", "structures"):
        raise click.UsageError("--af は structure=/structures= でのみ使用できます。")

    if data_type == "activities":
        target_chembl_id = resolve_chembl_target_id(value)
        records = fetch_activities(target_chembl_id)
        aggregated = activities.standardize_and_aggregate(records)
        activities.write_activities_tsv(aggregated, output)
    elif data_type == "structure":
        if af:
            fetch_af_structure(resolve_uniprot_accession(value), output, fmt=fmt)
        else:
            fetch_structure(value, output, fmt=fmt)
    elif data_type == "structures":
        if fmt is None:
            raise click.UsageError("structures=... では --type (cif/pdb) の指定が必須です。")
        identifiers = [x.strip() for x in value.split(",") if x.strip()]
        if not identifiers:
            raise click.UsageError(f"構造IDが指定されていません: {spec!r}")
        if af:
            accessions = [resolve_uniprot_accession(identifier) for identifier in identifiers]
            fetch_af_structures(accessions, output, fmt.lower())
        else:
            fetch_structures(identifiers, output, fmt.lower())
    else:
        raise click.UsageError(
            f"未対応のデータ種別です: {data_type!r} (対応: {', '.join(SUPPORTED_TYPES)})"
        )

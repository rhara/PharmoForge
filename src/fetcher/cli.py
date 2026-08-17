from pathlib import Path

import click

from afdb import fetch_structure as fetch_af_structure
from chembl import fetch_activities
from core.logging_utils import get_logger
from idmap import (
    looks_like_uniprot_accession,
    pdb_id_to_uniprot_accessions,
    resolve_chembl_target_id,
    resolve_uniprot_accession,
)
from rcsb import fetch_structure

from . import activities

logger = get_logger(__name__)

SUPPORTED_TYPES = ("activity", "structure")


@click.command("fetch")
@click.argument("spec")
@click.option(
    "--af",
    is_flag=True,
    default=False,
    help="Fetch structure= from AlphaFold DB instead of the default RCSB PDB. "
    "The identifier must be a UniProt entry name/accession.",
)
@click.option(
    "--outdir",
    "-o",
    "outdir",
    required=True,
    type=click.Path(path_type=Path, file_okay=False),
    help="Output directory. The file name is derived automatically from the identifier.",
)
@click.option(
    "--type",
    "fmt",
    type=click.Choice(["cif", "pdb", "fasta"], case_sensitive=False),
    default="cif",
    show_default=True,
    help="Structure file format (for structure=).",
)
def fetch_cmd(spec: str, af: bool, outdir: Path, fmt: str):
    """Fetch data. SPEC has the form <data type>=<identifier>.

    The output always goes to the directory given by --outdir/-o, with the file name derived
    automatically from the identifier (activity= produces `<identifier>_activity.tsv`,
    structure= produces `<identifier>.<fmt>`. With --af (except --type=fasta), the file name
    uses the UniProt accession with `_AF` appended).

    --type=fasta always fetches the sequence directly from UniProt itself and saves it as
    `<UniProt accession>.fasta` (no `_AF` suffix). --af is not needed: if the identifier is a
    PDB ID (e.g. `9CSK`), the UniProt accession(s) linked to that PDB entry are resolved
    automatically (multiple for complexes); a UniProt entry name/accession (e.g. `R1AB_SARS2`,
    `P61626`) is resolved as-is.

    \b
    Examples:
      pf fetch activity=CDK4_HUMAN --outdir data
      pf fetch structure=9CSK --type=cif --outdir data
      pf fetch structure=6P8F,7SJ3,9CSK --type pdb --outdir data
      pf fetch structure=9CSK --type=fasta --outdir data
      pf fetch structure=R1AB_SARS2 --type=fasta --outdir data
      pf fetch structure=P61626 --af --type=cif --outdir data
      pf fetch structure=TYK2_HUMAN --af --type=cif --outdir data
      pf fetch structure=P61626,CDK4_HUMAN --af --type pdb --outdir data
    """
    if "=" not in spec:
        raise click.UsageError(
            f"SPECは <データ種別>=<識別子> の形式で指定してください(例: activity=CDK4_HUMAN): {spec!r}"
        )
    data_type, _, value = spec.partition("=")
    data_type = data_type.strip().lower()
    value = value.strip()

    if af and data_type != "structure":
        raise click.UsageError("--af は structure= でのみ使用できます。")

    if data_type == "activity":
        target_chembl_id = resolve_chembl_target_id(value)
        records = fetch_activities(target_chembl_id)
        aggregated = activities.standardize_and_aggregate(records)
        activities.write_activities_tsv(aggregated, outdir / f"{value}_activity.tsv")
    elif data_type == "structure":
        identifiers = [x.strip() for x in value.split(",") if x.strip()]
        if not identifiers:
            raise click.UsageError(f"構造IDが指定されていません: {spec!r}")
        if fmt == "fasta":
            if af:
                logger.warning(
                    "--af has no effect on FASTA output (fasta is always fetched directly "
                    "from UniProt, and --af is not required to resolve entry names either)."
                )
            accessions: list[str] = []
            for identifier in identifiers:
                if af or looks_like_uniprot_accession(identifier) or "_" in identifier:
                    candidates = [resolve_uniprot_accession(identifier)]
                else:
                    candidates = pdb_id_to_uniprot_accessions(identifier)
                for accession in candidates:
                    if accession not in accessions:
                        accessions.append(accession)
            for accession in accessions:
                fetch_af_structure(accession, outdir / f"{accession}.fasta", fmt="fasta")
        else:
            for identifier in identifiers:
                if af:
                    accession = resolve_uniprot_accession(identifier)
                    fetch_af_structure(accession, outdir / f"{accession}_AF.{fmt}", fmt=fmt)
                else:
                    fetch_structure(identifier, outdir / f"{identifier.upper()}.{fmt}", fmt=fmt)
    else:
        raise click.UsageError(
            f"未対応のデータ種別です: {data_type!r} (対応: {', '.join(SUPPORTED_TYPES)})"
        )

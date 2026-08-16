from pathlib import Path

import click

from afdb import fetch_structure as fetch_af_structure
from chembl import fetch_activities
from idmap import resolve_chembl_target_id, resolve_uniprot_accession
from rcsb import fetch_structure

from . import activities

SUPPORTED_TYPES = ("activity", "structure")


@click.command("fetch")
@click.argument("spec")
@click.option(
    "--af",
    is_flag=True,
    default=False,
    help="structure=の取得元をAlphaFold DBにする(既定はRCSB PDB)。"
    "識別子はUniProt entry name/accessionを指定する。",
)
@click.option(
    "--outdir",
    "-o",
    "outdir",
    required=True,
    type=click.Path(path_type=Path, file_okay=False),
    help="出力先ディレクトリ。ファイル名は識別子から自動的に決まる。",
)
@click.option(
    "--type",
    "fmt",
    type=click.Choice(["cif", "pdb"], case_sensitive=False),
    default="cif",
    show_default=True,
    help="構造ファイルのフォーマット(structure=用)。",
)
def fetch_cmd(spec: str, af: bool, outdir: Path, fmt: str):
    """データを取得する。SPECは <データ種別>=<識別子> の形式。

    出力先は常に--outdir/-oで指定するディレクトリで、ファイル名は識別子から自動的に決まる
    (activity=は`<識別子>_activity.tsv`、structure=は`<識別子>.<fmt>`。--af指定時はUniProt
    accessionに`_AF`を付けたものをファイル名に使う)。

    \b
    例:
      pf fetch activity=CDK4_HUMAN --outdir data
      pf fetch structure=9CSK --type=cif --outdir data
      pf fetch structure=6P8F,7SJ3,9CSK --type pdb --outdir data
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

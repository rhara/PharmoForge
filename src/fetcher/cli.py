from pathlib import Path

import click

from chembl import fetch_activities
from idmap import resolve_chembl_target_id
from rcsb import fetch_structure, fetch_structures

from . import activities

SUPPORTED_TYPES = ("activities", "structure", "structures")


@click.command("fetch")
@click.argument("spec")
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
def fetch_cmd(spec: str, output: Path, fmt: str | None):
    """データを取得する。SPECは <データ種別>=<識別子> の形式。

    \b
    例:
      pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
      pf fetch structure=9CSK --type=cif --output data/9csk.cif
      pf fetch structures=6P8F,7SJ3,9CSK --type pdb --output data
    """
    if "=" not in spec:
        raise click.UsageError(
            f"SPECは <データ種別>=<識別子> の形式で指定してください(例: activities=CDK4_HUMAN): {spec!r}"
        )
    data_type, _, value = spec.partition("=")
    data_type = data_type.strip().lower()
    value = value.strip()

    if data_type == "activities":
        target_chembl_id = resolve_chembl_target_id(value)
        records = fetch_activities(target_chembl_id)
        aggregated = activities.standardize_and_aggregate(records)
        activities.write_activities_tsv(aggregated, output)
    elif data_type == "structure":
        fetch_structure(value, output, fmt=fmt)
    elif data_type == "structures":
        if fmt is None:
            raise click.UsageError("structures=... では --type (cif/pdb) の指定が必須です。")
        pdb_ids = [x.strip() for x in value.split(",") if x.strip()]
        if not pdb_ids:
            raise click.UsageError(f"構造IDが指定されていません: {spec!r}")
        fetch_structures(pdb_ids, output, fmt.lower())
    else:
        raise click.UsageError(
            f"未対応のデータ種別です: {data_type!r} (対応: {', '.join(SUPPORTED_TYPES)})"
        )

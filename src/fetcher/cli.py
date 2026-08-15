from pathlib import Path

import click

from idmap import resolve_chembl_target_id

from . import chembl, structures

SUPPORTED_TYPES = ("activities", "structure")


@click.command("fetch")
@click.argument("spec")
@click.option(
    "--output", "-o", "output", required=True, type=click.Path(path_type=Path), help="出力ファイルパス"
)
def fetch_cmd(spec: str, output: Path):
    """データを取得する。SPECは <データ種別>=<識別子> の形式。

    \b
    例:
      pf fetch activities=CDK4_HUMAN --output data/cdk4_human_activities.tsv
      pf fetch structure=9CSK --output data/9csk.cif
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
        records = chembl.fetch_activities(target_chembl_id)
        chembl.write_activities_tsv(records, output)
    elif data_type == "structure":
        structures.fetch_structure(value, output)
    else:
        raise click.UsageError(
            f"未対応のデータ種別です: {data_type!r} (対応: {', '.join(SUPPORTED_TYPES)})"
        )

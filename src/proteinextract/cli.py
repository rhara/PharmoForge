from pathlib import Path

import click

from .extract import extract_structure


@click.command("protein-extract")
@click.argument("input_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--chains",
    default=None,
    help="抽出するチェーンIDをカンマ区切りで指定(例: A,G,H,I,J)。省略時は全チェーン。",
)
@click.option(
    "--remove-water",
    is_flag=True,
    default=False,
    help="水分子(HOH等)を除去する",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    required=True,
    type=click.Path(path_type=Path),
    help="出力構造ファイルパス(拡張子.pdb/.cifで形式を判別)",
)
def protein_extract_cmd(input_path: Path, chains: str | None, remove_water: bool, output_path: Path):
    """構造ファイル(PDB/CIF)から指定チェーンを抽出し、必要に応じて水分子を除去する。

    入力・出力とも拡張子(`.pdb`/`.cif`)で形式を自動判別するため、異なる形式間でも変換できる。

    \b
    例:
      pf protein-extract data/cdk2/2CCH.cif --chains=A,G,H,I,J --remove-water --output data/cdk2/2CCH_main.cif
      pf protein-extract data/tyk2/6NZP.cif --remove-water --output data/tyk2/6NZP_nowater.pdb
    """
    chain_list = [c.strip() for c in chains.split(",") if c.strip()] if chains else None
    extract_structure(input_path, output_path, chains=chain_list, remove_water=remove_water)

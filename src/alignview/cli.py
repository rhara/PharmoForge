from pathlib import Path

import click

from .view import launch_alignment_view


@click.command("align-view")
@click.argument("structure_paths", nargs=-1, required=True, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--method",
    type=click.Choice(["align", "super", "cealign", "number"]),
    default="align",
    show_default=True,
    help="アラインメント手法(align/number=配列非依存で残基番号ベース、super/cealign=配列非依存で構造ベース)",
)
@click.option(
    "--align-margin",
    type=int,
    default=20,
    show_default=True,
    help="--method align時、誤対応防止のためtargetの探索範囲をmobileの残基番号レンジ±この値に絞り込む",
)
@click.option(
    "--pymol-env",
    default="pymol",
    show_default=True,
    help="PyMOLをインストールした専用conda/mamba環境名",
)
def align_view_cmd(structure_paths: tuple[Path, ...], method: str, align_margin: int, pymol_env: str):
    """複数のPDB/CIF構造をPyMOLで開き、先頭の構造に他を重ね合わせる(アラインメント)。

    同一蛋白の構造間ではPDBの残基番号が(UniProt基準等で)揃っている前提のもと、
    `align`/`number`はいずれも配列アラインメントではなく残基番号の対応付けを利用する。

    \b
    - align (既定): targetの探索範囲をmobileの残基番号レンジ±`--align-margin`に絞り込んだ上で
      PyMOLの配列ベースalignを行う(絞り込みで対応が取れない場合は絞り込みなしにフォールバック)。
      外れ値を反復的に除外して精密化するため、対応する範囲では最も高精度になりやすい。
    - number: 配列アラインメントを一切行わず、残基番号が一致するCA原子同士を直接対応付けて
      剛体重ね合わせを行う(structfit/ProDyによる計算)。最も直接的・決定的な対応付け。
    - super / cealign: 配列非依存の構造ベースアラインメント。異なる蛋白同士の比較や、構造間で
      残基番号が揃っていない場合に使う。

    \b
    例:
      pf align-view data/TYK2_HUMAN_af.cif data/6NZP.cif data/4OLI.cif data/5C03.cif
      pf align-view data/TYK2_HUMAN_af.cif data/6NZP.cif --method number
      pf align-view data/9csk.cif data/1abc.cif --method cealign
    """
    launch_alignment_view(list(structure_paths), method=method, pymol_env=pymol_env, align_margin=align_margin)

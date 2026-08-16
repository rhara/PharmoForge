from pathlib import Path

import click

from .view import launch_alignment_view

# 拡張子省略時に試す順(優先度順)
_AUTO_EXTENSIONS = (".cif", ".pdb")


def _resolve_structure_paths(tokens: tuple[str, ...]) -> list[Path]:
    """`--indir`と拡張子省略に対応しつつ、tokens列を構造ファイルパスのリストに解決する。"""
    paths: list[Path] = []
    indir: Path | None = None
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "--indir":
            i += 1
            if i >= len(tokens):
                raise click.UsageError("--indir にはディレクトリを指定してください")
            indir = Path(tokens[i])
            i += 1
            continue
        paths.append(_resolve_one(token, indir))
        i += 1
    if not paths:
        raise click.UsageError("構造ファイルを1つ以上指定してください")
    return paths


def _resolve_one(token: str, indir: Path | None) -> Path:
    raw = Path(token)
    # "/"を含む(または絶対パスの)指定は--indirによらずそのまま(カレント相対 or 絶対)扱う
    if indir is None or "/" in token or raw.is_absolute():
        base = raw
    else:
        base = indir / raw

    if base.suffix:
        if not base.exists():
            raise click.UsageError(f"ファイルが見つかりません: {base}")
        return base

    for ext in _AUTO_EXTENSIONS:
        candidate = base.with_suffix(ext)
        if candidate.exists():
            return candidate
    tried = ", ".join(str(base.with_suffix(ext)) for ext in _AUTO_EXTENSIONS)
    raise click.UsageError(f"ファイルが見つかりません({tried})")


@click.command("align-view", context_settings={"ignore_unknown_options": True})
@click.argument("tokens", nargs=-1, type=click.UNPROCESSED)
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
def align_view_cmd(tokens: tuple[str, ...], method: str, align_margin: int, pymol_env: str):
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
    --indir DIR: 繰り返し指定可能。以降のファイル名(拡張子省略可、.cif優先、次に.pdb)を
      DIR配下から解決する。"/"を含む指定(または絶対パス)は--indirによらず
      カレントディレクトリ相対 or 絶対パスとして扱う。

    \b
    例:
      pf align-view data/TYK2_HUMAN_af.cif data/6NZP.cif data/4OLI.cif data/5C03.cif
      pf align-view --indir data/cyp P08604_AF 1PQ2_ad 3IBD_abcde 3NXU_abh
      pf align-view --indir data/cyp P08604_AF 1PQ2_ad --indir data/other 9XYZ
      pf align-view data/TYK2_HUMAN_af.cif data/6NZP.cif --method number
    """
    structure_paths = _resolve_structure_paths(tokens)
    launch_alignment_view(structure_paths, method=method, pymol_env=pymol_env, align_margin=align_margin)

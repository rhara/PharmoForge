from pathlib import Path

import click

from structio.resolve import resolve_structure_tokens

from .report import build_report, load_labeled_structures


@click.command("sequence-align", context_settings={"ignore_unknown_options": True})
@click.argument("tokens", nargs=-1, type=click.UNPROCESSED)
@click.option(
    "--reference",
    default=None,
    help=(
        "残基置換一覧の基準。'ラベル:チェーンID'(例: P24941_AF:A)で読み込んだ構造の"
        "1チェーンを指定するか、アミノ酸配列(1文字表記)を直接指定できる。省略時は置換一覧を出力しない。"
    ),
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="レポートの保存先(省略時は標準出力)",
)
def sequence_align_cmd(tokens: tuple[str, ...], reference: str | None, output_path: Path | None):
    """複数のPDB/CIF構造から蛋白配列を抽出し、FASTA・pairwise identity・
    (--reference指定時)基準チェーンに対する残基置換一覧を出力する。

    配列はチェーンごとにCA原子(観測された残基のみ)から抽出するため、電子密度が
    見えず欠損した残基は含まれない(UniProtの完全配列とは異なりうる)。

    --referenceで基準を指定した場合、残基置換一覧を出力する。基準の指定方法は2通り:

    \b
    - 'ラベル:チェーンID'(例: P24941_AF:A): 読み込んだ構造の1チェーンを基準にする。
      残基番号ベースの対応付けのみを用いる(同一蛋白の構造間ではPDBの残基番号が
      揃っている前提。`pf align-view --method number`と同じ前提)。番号体系が
      異なる構造間では「対応が取れませんでした」と表示される。
    - アミノ酸配列(1文字表記、コロンを含まない): 構造を伴わない任意配列
      (UniProt正規配列やユーザー指定の基準配列)を基準にする。この場合は配列
      アラインメントを用いるため、残基番号が揃っていない構造間でも比較できる。

    \b
    --indir DIR: 繰り返し指定可能。以降のファイル名(拡張子省略可、.cif優先、
      次に.pdb)をDIR配下から解決する。"/"を含む指定(または絶対パス)は
      --indirによらずカレントディレクトリ相対 or 絶対パスとして扱う(align-viewと同様)。

    \b
    例:
      pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
      pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a --reference P24941_AF:A
      pf sequence-align data/braf/P15056_AF.cif data/braf/3OG7_ac.cif --reference P15056_AF:A -o report.txt
      pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab 1HCL_a
    """
    structure_paths = resolve_structure_tokens(tokens)
    structures = load_labeled_structures(structure_paths)
    report = build_report(structures, reference)
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report)
        click.echo(f"saved report to {output_path}")
    else:
        click.echo(report)

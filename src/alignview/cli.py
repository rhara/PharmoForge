import click

from structio.resolve import resolve_structure_tokens

from .view import launch_alignment_view


@click.command("align-view", context_settings={"ignore_unknown_options": True})
@click.argument("tokens", nargs=-1, type=click.UNPROCESSED)
@click.option(
    "--method",
    type=click.Choice(["align", "super", "cealign", "number"]),
    default="align",
    show_default=True,
    help="Alignment method (align/number = sequence-independent, residue-number based; super/cealign = sequence-independent, structure based)",
)
@click.option(
    "--align-margin",
    type=int,
    default=20,
    show_default=True,
    help="With --method align, restrict target's search range to mobile's residue-number range +/- this margin, to avoid mismatches",
)
@click.option(
    "--pymol-env",
    default="pymol",
    show_default=True,
    help="Name of the dedicated conda/mamba environment PyMOL is installed in",
)
def align_view_cmd(tokens: tuple[str, ...], method: str, align_margin: int, pymol_env: str):
    """Open multiple PDB/CIF structures in PyMOL and superpose the rest onto the first (alignment).

    Assuming PDB residue numbers are consistent (e.g. UniProt-based) across structures of the
    same protein, both `align` and `number` use residue-number correspondence rather than
    sequence alignment.

    \b
    - align (default): restricts target's search range to mobile's residue-number range +/-
      `--align-margin`, then runs PyMOL's sequence-based align (falls back to unrestricted align
      if the restricted range yields no matches). Iteratively excludes outliers to refine the fit,
      so it tends to be the most accurate within the matched range.
    - number: skips sequence alignment entirely and directly pairs CA atoms with matching residue
      numbers for a rigid-body superposition (computed via structfit/ProDy). The most direct,
      deterministic correspondence.
    - super / cealign: sequence-independent, structure-based alignment. Use for comparing
      different proteins, or when residue numbers aren't consistent across structures.

    \b
    --indir DIR: can be repeated. Subsequent file names (extension optional, .cif preferred,
      then .pdb) are resolved under DIR. A token containing "/" (or an absolute path) is treated
      as current-directory-relative or absolute regardless of --indir.

    \b
    Examples:
      pf align-view data/TYK2_HUMAN_af.cif data/6NZP.cif data/4OLI.cif data/5C03.cif
      pf align-view --indir data/cyp P08604_AF 1PQ2_ad 3IBD_abcde 3NXU_abh
      pf align-view --indir data/cyp P08604_AF 1PQ2_ad --indir data/other 9XYZ
      pf align-view data/TYK2_HUMAN_af.cif data/6NZP.cif --method number
    """
    structure_paths = resolve_structure_tokens(tokens)
    launch_alignment_view(structure_paths, method=method, pymol_env=pymol_env, align_margin=align_margin)

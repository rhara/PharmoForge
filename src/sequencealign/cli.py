from pathlib import Path

import click

from structio.resolve import resolve_structure_tokens

from .report import DEFAULT_ALIGN_WIDTH, build_report, load_labeled_structures


@click.command("sequence-align", context_settings={"ignore_unknown_options": True})
@click.argument("tokens", nargs=-1, type=click.UNPROCESSED)
@click.option(
    "--reference",
    default=None,
    help=(
        "Reference for the substitution list: either 'label:chain_id' (e.g. P24941_AF:A) "
        "identifying a chain from a loaded structure, or an amino acid sequence (one-letter code) "
        "given directly. Omit to skip the substitution list."
    ),
)
@click.option(
    "--width",
    type=int,
    default=DEFAULT_ALIGN_WIDTH,
    show_default=True,
    help="Number of residues per line in the residue-number alignment section.",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Path to save the report to (default: print to stdout).",
)
def sequence_align_cmd(tokens: tuple[str, ...], reference: str | None, width: int, output_path: Path | None):
    """Extract protein sequences from multiple PDB/CIF structures and report FASTA,
    pairwise identity, a residue-number alignment, and (with --reference) substitutions.

    Sequences are extracted per chain from observed CA atoms only, so residues not
    resolved in the electron density are excluded (this can differ from the full
    UniProt sequence).

    With --reference given, a substitution list is also printed. Two ways to specify it:

    \b
    - 'label:chain_id' (e.g. P24941_AF:A): use one chain from a loaded structure as the
      reference. Only residue-number-based correspondence is used (assumes PDB residue
      numbers line up across structures, same assumption as `pf align-view --method number`).
      Structures with a different numbering scheme are reported as "no correspondence found".
    - An amino acid sequence (one-letter code, no colon): use an arbitrary sequence not tied
      to any structure (e.g. a UniProt canonical sequence) as the reference. This uses a
      sequence alignment, so it also works across structures whose residue numbers don't align.

    \b
    --indir DIR: repeatable. Resolves following filenames (extension optional, .cif
      tried first, then .pdb) under DIR. A token containing "/" (or an absolute path)
      is used as-is (relative to cwd or absolute), regardless of --indir (same as align-view).

    \b
    Examples:
      pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
      pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a --reference P24941_AF:A
      pf sequence-align data/braf/P15056_AF.cif data/braf/3OG7_ac.cif --reference P15056_AF:A -o report.txt
      pf sequence-align --reference MENFQKV...PHLRL --indir data/cdk2 1AQ1_ab 1HCL_a
      pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160
    """
    structure_paths = resolve_structure_tokens(tokens)
    structures = load_labeled_structures(structure_paths)
    report = build_report(structures, reference, align_width=width)
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report)
        click.echo(f"saved report to {output_path}")
    else:
        click.echo(report)

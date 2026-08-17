from pathlib import Path

import click

from structio.resolve import resolve_structure_tokens

from .report import DEFAULT_ALIGN_WIDTH, build_report, load_labeled_structures

# .fasta も入力として許容する(sequence-alignは配列比較のみで3次元構造を必要としないため)。
_INPUT_EXTENSIONS = (".cif", ".mmcif", ".pdb", ".fasta")


@click.command("sequence-align", context_settings={"ignore_unknown_options": True})
@click.argument("tokens", nargs=-1, type=click.UNPROCESSED)
@click.option(
    "--method",
    type=click.Choice(["number", "align"]),
    default="align",
    show_default=True,
    help=(
        "Alignment method for the alignment section. 'align' (default) uses pairwise "
        "sequence alignment (Biopython PairwiseAligner) against the first input's first "
        "chain, so it works even when residue numbering differs across inputs (e.g. a "
        "full-length sequence vs. a domain-only structure). 'number' instead lines "
        "sequences up by raw PDB residue number, which is only meaningful when numbering "
        "is already consistent across inputs (e.g. same UniProt numbering)."
    ),
)
@click.option(
    "--width",
    type=int,
    default=DEFAULT_ALIGN_WIDTH,
    show_default=True,
    help="Number of residues per line in the alignment section.",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Path to save the report to (default: print to stdout).",
)
def sequence_align_cmd(tokens: tuple[str, ...], method: str, width: int, output_path: Path | None):
    """Extract protein sequences from multiple PDB/CIF structures (and/or plain FASTA
    files) and report pairwise identity and an alignment.

    Sequences from PDB/CIF are extracted per chain from observed CA atoms only, so
    residues not resolved in the electron density are excluded (this can differ from
    the full UniProt sequence). A .fasta input has no atomic structure, so its
    sequence(s) are used as-is and it is excluded from pairwise identity (which needs
    both sides' atoms). To add a reference sequence (e.g. a UniProt canonical
    sequence) to the alignment, just include it as a regular .fasta input token
    alongside the structures.

    --method controls how the alignment section lines sequences up (see --method
    above): 'align' (default) uses real sequence alignment, so numbering doesn't need
    to correspond across inputs (e.g. --indir data/mpro P0DTD1 6LU7_abc: the
    full-length polyprotein FASTA numbered 1..7096 vs. the Mpro domain structure
    numbered locally 1..306 — 'number' would wrongly line these up by raw position
    instead of by homology). Use --method number only when residue numbering is
    already known to be consistent across inputs (e.g. the same UniProt numbering).

    \b
    --indir DIR: repeatable. Resolves following filenames (extension optional, tried
      in order .cif, .mmcif, .pdb, .fasta) under DIR. A token containing "/" (or an
      absolute path) is used as-is (relative to cwd or absolute), regardless of
      --indir (same as align-view).

    \b
    Examples:
      pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a
      pf sequence-align --indir data/braf P15056_AF 4MNF_ac --width 160 -o report.txt
      pf sequence-align --indir data/mpro P0DTD1.fasta 6LU7_abc
      pf sequence-align --indir data/mpro P0DTD1 6LU7_abc
      pf sequence-align --indir data/cdk2 P24941_AF 1AQ1_ab 1HCL_a --method number
    """
    structure_paths = resolve_structure_tokens(tokens, extensions=_INPUT_EXTENSIONS)
    structures = load_labeled_structures(structure_paths)
    report = build_report(structures, align_width=width, method=method)
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report)
        click.echo(f"saved report to {output_path}")
    else:
        click.echo(report)

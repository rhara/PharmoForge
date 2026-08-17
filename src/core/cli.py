import click

from alignview.cli import align_view_cmd
from fetcher.cli import fetch_cmd
from pocketfinder.cli import find_pocket_cmd, show_pocket_cmd, view_pocket_cmd
from proteinanalyzer.cli import protein_info_cmd
from proteinextract.cli import protein_extract_cmd
from proteinprep.cli import prep_protein_cmd
from scaffoldanalyzer.cli import analyze_scaffolds_cmd
from sequencealign.cli import sequence_align_cmd


@click.group()
def cli():
    """pf: PharmoForge unified CLI."""


cli.add_command(fetch_cmd)
cli.add_command(analyze_scaffolds_cmd)
cli.add_command(prep_protein_cmd)
cli.add_command(protein_info_cmd)
cli.add_command(align_view_cmd)
cli.add_command(protein_extract_cmd)
cli.add_command(sequence_align_cmd)
cli.add_command(find_pocket_cmd)
cli.add_command(view_pocket_cmd)
cli.add_command(show_pocket_cmd)


if __name__ == "__main__":
    cli()

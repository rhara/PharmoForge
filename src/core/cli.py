import click

from fetcher.cli import fetch_cmd
from proteinanalyzer.cli import protein_info_cmd
from proteinprep.cli import prep_protein_cmd
from scaffoldanalyzer.cli import analyze_scaffolds_cmd


@click.group()
def cli():
    """pf: PharmoForge unified CLI."""


cli.add_command(fetch_cmd)
cli.add_command(analyze_scaffolds_cmd)
cli.add_command(prep_protein_cmd)
cli.add_command(protein_info_cmd)


if __name__ == "__main__":
    cli()

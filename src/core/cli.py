import click

from fetcher.cli import fetch_cmd


@click.group()
def cli():
    """pf: PharmoForge unified CLI."""


cli.add_command(fetch_cmd)


if __name__ == "__main__":
    cli()

"""biome-finder command line."""

import click
from . import survey as survey_mod


@click.group()
def main() -> None:
    """Find and organize public INSDC data for a target biome."""


@main.command()
@click.option(
    "--biomes",
    "biome_paths",
    multiple=True,
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Biome YAML config (repeatable).",
)
@click.argument("extra", nargs=-1, metavar="[BIOME_YAML]...", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--out",
    type=click.Path(file_okay=False),
    default="survey",
    show_default=True,
    help="Output folder: survey.tsv, biomes/<name>/*.tsv, report.html.",
)
def survey(biome_paths: tuple[str, ...], extra: tuple[str, ...], out: str) -> None:
    """Count public INSDC data per candidate biome.

    Configs can also be given as extra arguments, so `--biomes biomes/*.yaml` works."""
    try:
        survey_mod.run(biome_paths + extra, out)
    except ValueError as e:  # invalid biome config
        raise click.ClickException(str(e)) from e


if __name__ == "__main__":
    main()

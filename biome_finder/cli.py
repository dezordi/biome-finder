"""biome-finder command line."""

import click
from . import log
from . import search as search_mod
from . import survey as survey_mod


@click.group()
@click.option("-v", "--verbose", is_flag=True, help="Also log every HTTP request (cached or not).")
def main(verbose: bool) -> None:
    """Find and organize public INSDC data for a target biome."""
    log.setup(verbose)


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


@main.command()
@click.option("--biome", "biome_path", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--survey",
    "survey_dir",
    type=click.Path(file_okay=False),
    default="survey",
    show_default=True,
    help="Survey folder holding biomes/<name>/records.tsv.",
)
@click.option("--out", type=click.Path(dir_okay=False), default="samples.tsv", show_default=True)
@click.option("--since", default="", metavar="YYYY-MM-DD", help="Only records public on or after this date.")
@click.option(
    "--taxdump", default="", metavar="YYYY-MM-DD", help="NCBI taxonomy release of that month (default: the latest)."
)
@click.option(
    "--gold",
    "gold_export",
    default="",
    type=click.Path(dir_okay=False),
    help="GOLD data export (goldData.xlsx) for GOLD ecosystem evidence; see the search docs to download it.",
)
def search(biome_path: str, survey_dir: str, out: str, since: str, taxdump: str, gold_export: str) -> None:
    """Build samples.tsv for one biome from its survey records (merged into an existing one)."""
    try:
        search_mod.run(biome_path, survey_dir, out, since, taxdump, gold_export)
    except ValueError as e:  # survey not run for this biome, or no such taxdump release
        raise click.ClickException(str(e)) from e


if __name__ == "__main__":
    main()

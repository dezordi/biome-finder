"""Survey dashboard: one self-contained `report.html`, built with Panel and plotly.

Drawn only from the tables in a survey folder, so it can be redrawn without querying ENA.
Two tabs: an overview of all catalogues, and one catalogue picked from a dropdown. Every
catalogue view is precomputed into the page (`save(embed=True)`), so the file works offline.

Everything shown is for the catalogue's kingdom: its genomes, and only the assemblies and runs
with evidence for that kingdom (`kingdom_match` in records.tsv). The overview compares many
catalogues, so it uses heatmaps with the count in each cell; a catalogue view uses bars.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import panel as pn
import plotly.express as px
import plotly.graph_objects as go
from bokeh.resources import INLINE
from plotly.subplots import make_subplots

from .constants import LONG_READ_PLATFORMS, MIN_READS

SUBMITTED, BY_MGNIFY = "MAGs (submitters)", "MAGs (MGnify)"
ASSEMBLIES, RUNS = "Assemblies", "Runs ≥1M"
KINDS = [SUBMITTED, BY_MGNIFY, ASSEMBLIES, RUNS]
ROUTES = {
    "MAGs to add": "→ genomes-catalogue-pipeline",
    "Samples to bin": "→ genomes-generation",
    "Samples to assemble": "→ miassembler, then genomes-generation",
}
TOP_STUDIES = 15
PROJECTS = (
    "Projects are ENA/INSDC BioProjects (PRJEB, PRJNA, PRJDB): one submission grouping the samples "
    "and data of one research effort, usually one paper. Counted per route: the projects holding "
    "the MAGs, the MAGs to add, and the assemblies or runs of the samples to bin or assemble. "
    "Many projects mean more varied data, and "
    "more to curate, since curation reads each project's paper and metadata. For MAGs (MGnify), "
    "the projects are MGnify's own analysis projects, not the original ones."
)
ACCENT = "#0072B2"
COLORS = {
    SUBMITTED: "#0072B2",
    BY_MGNIFY: "#56B4E9",
    ASSEMBLIES: "#E69F00",
    RUNS: "#009E73",
    "MAGs": "#0072B2",
    "MAGs to add": "#0072B2",
    "Assemblies to bin": "#E69F00",
    "Runs to assemble": "#009E73",
    "Samples to bin": "#E69F00",
    "Samples to assemble": "#009E73",
    "short reads": "#009E73",
    "long reads": "#CC79A7",
    "in catalogue": "#BBBBBB",
    "Bacteria": "#0072B2",
    "Archaea": "#D55E00",
    "Eukaryota": "#009E73",
}
# column: (title, tooltip)
COLUMNS = {
    "biome": (
        "Catalogue",
        "The biome config: one catalogue per biome and kingdom, as in MGnify. *: MGnify already has "
        "this catalogue, so the survey counts an update",
    ),
    "metagenome_taxon": ("Metagenome taxon", "ENA metagenome taxon the records were found under"),
    "tax_id": ("Tax ID", "NCBI tax ID of the metagenome taxon"),
    "kingdom": ("Kingdom", "prokaryotes (Bacteria + Archaea) or eukaryotes"),
    "mags": (
        "MAGs",
        "MAGs of the catalogue's kingdom in ENA: deposited by their authors, or produced by MGnify from "
        "this biome's public data (ENA stores most of those as 'binned metagenome' records)",
    ),
    "to_add": (
        "MAGs to add",
        "MAGs (both kinds) not in the existing MGnify catalogue (all of them, for a new "
        "catalogue): ready for genomes-catalogue-pipeline",
    ),
    "samples_bin": (
        "Samples to bin",
        "Samples with a metagenome assembly but no genome yet: input for genomes-generation. Only "
        "assemblies with evidence for the kingdom count (a genome of that kingdom from the same sample, "
        "or a project that names it). A sample can have several assemblies, e.g. by its submitter and "
        "by MGnify; hover the routes chart for the assembly count",
    ),
    "samples_assemble": (
        "Samples to assemble",
        "Samples with only reads (no genome, no assembly): input for miassembler, then "
        "genomes-generation. Only metagenomic WGS runs with at least 1 million reads and evidence for "
        "the kingdom count. A sample can have several runs (lanes, re-sequencing); hover the routes "
        "chart for the run count",
    ),
    "tbases": ("Runs Tb", "Terabases (10^12 bases, not bytes) in the runs to assemble"),
}
FIGURE_LAYOUT = {
    "template": "plotly_white",
    "height": 380,
    "margin": {"l": 10, "r": 10, "t": 36, "b": 10},
    "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0, "title": None},
    "font": {"size": 12},
}
TABLE_CSS = """
:host .tabulator-col-title { white-space: normal; word-break: normal; overflow-wrap: normal; text-overflow: clip; }
:host .tabulator-col { font-size: 12px; }
:host .tabulator-col-sorter, :host .tabulator-arrow { display: none; }  /* click the title to sort */
"""
ABOUT = """
Every number is for the catalogue's kingdom: its genomes, and the assemblies and runs with
evidence for that kingdom. Counts are proxies: a record belongs to a catalogue through the
metagenome taxon its submitter chose (plus the host, for host-associated biomes), not through
curation. Hover a column title for its definition.

**Where new genomes would come from.** Each sample enters at its most processed public data:
*MAGs to add* are ready for the catalogue; *samples to bin* have assemblies that need binning
first; *samples to assemble* have only reads, which need assembly and binning. Many MAGs to add
means the catalogue can grow now; many samples to bin or assemble means it would pay to generate
new MAGs. Samples, not records, because a sample can have several assemblies or runs.
"""


def projects(out: Path) -> pd.DataFrame:
    """Every catalogue's studies.tsv (records per route and project), with a `biome` column. Counts
    are made numeric because a catalogue with nothing on a route has an empty file, read as text."""
    studies = pd.concat(
        [
            pd.read_csv(p, sep="\t", dtype={"study_accession": str, "from_mgnify": str}).assign(biome=p.parent.name)
            for p in sorted((out / "biomes").glob("*/studies.tsv"))
        ],
        ignore_index=True,
    )
    counts = ["mags_to_add", "assemblies_to_bin", "samples_to_bin", "runs_to_assemble", "samples_to_assemble"]
    studies[counts] = studies[counts].apply(pd.to_numeric).fillna(0).astype(int)
    return studies


def load(out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """survey.tsv plus every catalogue's records.tsv, with a `biome` column and a `kind`
    (MAGs from submitters or MGnify, assemblies, runs ≥1M) for every record counted for the
    catalogue's kingdom; everything else is dropped."""
    survey = pd.read_csv(
        out / "survey.tsv", sep="\t", dtype={"catalogue_updated": str, "tax_id": str}, dtype_backend="numpy_nullable"
    )
    records = pd.concat(
        [
            pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False).assign(biome=p.parent.name)
            for p in sorted((out / "biomes").glob("*/records.tsv"))
        ],
        ignore_index=True,
    )
    records["first_public"] = pd.to_datetime(records["first_public"], errors="coerce")
    for col in ["read_count", "base_count", "lat", "lon"]:
        records[col] = pd.to_numeric(records[col], errors="coerce")
    genome = records["data_type"].isin(["mag", "bin"])
    by_mgnify = (records["data_type"] == "bin") | (records["from_mgnify"] == "True")
    matched = records["kingdom_match"] == "True"
    records["kind"] = np.select(
        [
            genome & by_mgnify,
            genome,
            (records["data_type"] == "contigs") & matched,
            (records["data_type"] == "reads") & matched & (records["read_count"] >= MIN_READS),
        ],
        [BY_MGNIFY, SUBMITTED, ASSEMBLIES, RUNS],
        default="",
    )
    return survey, records[records["kind"] != ""].copy()


def totals(rows: pd.DataFrame) -> pd.DataFrame:
    """survey.tsv rows as the dashboard's kingdom-only columns, with the samples of each route."""
    n = rows.fillna({c: 0 for c in rows.columns if pd.api.types.is_numeric_dtype(rows[c])})
    return pd.DataFrame(
        {
            "biome": rows["biome"],
            "metagenome_taxon": rows["metagenome_taxon"],
            "tax_id": rows["tax_id"],
            "kingdom": rows["kingdom"],
            "mags": n["mags"] + n["mgnify_bins"],
            "to_add": n["mags_to_add"] + n["mgnify_bins_to_add"],
            "to_bin": rows["assemblies_to_bin"],
            "to_assemble": rows["runs_to_assemble"],
            "samples_bin": rows["samples_to_bin"],
            "samples_assemble": rows["samples_to_assemble"],
            "tbases": rows["gbases_to_assemble"] / 1000,
        }
    ).reset_index(drop=True)


# --- figure helpers ---


def figure(
    frame: pd.DataFrame, plot: Callable[..., go.Figure], layout: dict[str, Any] | None = None, **kw: Any
) -> go.Figure:
    """A styled plotly express figure, or a 'no data' placeholder for an empty frame."""
    if frame.empty:
        fig = go.Figure()
        fig.add_annotation(text="no data", showarrow=False, font={"size": 16, "color": "grey"})
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
    else:
        fig = plot(frame, **kw)
    return fig.update_layout(FIGURE_LAYOUT | (layout or {}))


def heatmap(matrix: pd.DataFrame, colorscale: str = "Blues", hover: pd.DataFrame | None = None) -> go.Figure:
    """Counts as a heatmap: log-scaled color, the count written in each cell, height to fit.
    `hover`, of the same shape, replaces the count in the hover label."""
    matrix = matrix.astype(float)
    text = matrix.map(lambda v: "" if pd.isna(v) or v == 0 else f"{v:,.0f}")
    hover = text if hover is None else hover
    fig = go.Figure(
        go.Heatmap(
            z=np.log10(matrix.where(matrix > 0)),
            x=[str(c) for c in matrix.columns],
            y=list(matrix.index),
            text=text.values,
            texttemplate="%{text}",
            colorscale=colorscale,
            showscale=False,
            xgap=2,
            ygap=2,
            customdata=hover.values,
            hovertemplate="%{y} · %{x}: %{customdata}<extra></extra>",
        )
    )
    fig.update_yaxes(autorange="reversed", ticksuffix="  ")
    fig.update_xaxes(side="top", type="category", tickangle=0)
    return fig.update_layout(
        FIGURE_LAYOUT | {"height": 70 + 26 * len(matrix), "margin": {"l": 10, "r": 10, "t": 50, "b": 10}}
    )


def depth_figure(records: pd.DataFrame, by: str) -> go.Figure:
    """Reads per run as boxes from precomputed percentiles (5/25/50/75/95)."""
    runs = records[records["kind"] == RUNS]
    runs = runs.assign(
        platform=runs["instrument_platform"].isin(LONG_READ_PLATFORMS).map({True: "long reads", False: "short reads"})
    )
    q = runs.groupby(["platform", by])["read_count"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
    fig = go.Figure(
        [
            go.Box(
                name=platform,
                y=list(g.index.get_level_values(by)),
                lowerfence=g[0.05],
                q1=g[0.25],
                median=g[0.5],
                q3=g[0.75],
                upperfence=g[0.95],
                orientation="h",
                marker_color=COLORS[platform],
            )
            for platform, g in q.groupby(level="platform")
        ]
    )
    n = max(runs[by].nunique(), 1)
    return (
        fig.update_layout(FIGURE_LAYOUT | {"height": max(260, 80 + 34 * n)}, boxmode="group")
        .update_xaxes(type="log", title="reads per run")
        .update_yaxes(categoryorder="category descending")
    )


def map_figure(records: pd.DataFrame, color: str, colors: dict[str, str]) -> go.Figure:
    """Assemblies (circles) and runs (diamonds) on a 1° grid; MAGs have no coordinates in ENA.
    Colored by `color`; when that isn't the kind, the legend lists each `color` value once and
    a click or double-click acts on all its markers."""
    keys = list(dict.fromkeys([color, "kind", "lat", "lon"]))
    data = (
        records[records["kind"].isin([ASSEMBLIES, RUNS])]
        .dropna(subset=["lat", "lon"])
        .assign(lat=lambda d: d["lat"].round(), lon=lambda d: d["lon"].round())
        .groupby(keys, as_index=False)
        .agg(records=("accession", "size"), projects=("study_accession", "nunique"))
    )
    fig = figure(
        data,
        px.scatter_geo,
        layout={"height": 520},
        lat="lat",
        lon="lon",
        size="records",
        color=color,
        symbol="kind",
        symbol_map={ASSEMBLIES: "circle", RUNS: "diamond"},
        color_discrete_map=colors,
        hover_data=["kind", "records", "projects"],
        size_max=22,
        opacity=0.65,
        projection="natural earth",
    )
    if color != "kind":  # one plain legend entry per group, whatever shapes it has
        for trace in fig.data:
            trace.update(legendgroup=trace.name.split(", ")[0], showlegend=False)
        for group in dict.fromkeys(t.legendgroup for t in fig.data):
            fig.add_scattergeo(
                lat=[None], lon=[None], mode="markers", name=group, legendgroup=group, marker_color=colors[group]
            )
        fig.update_layout(legend_groupclick="togglegroup")
    return kind_buttons(fig).update_geos(
        showcountries=True, countrycolor="#d0d0d0", showland=True, landcolor="#f4f4f4", showframe=False
    )


def kind_buttons(fig: go.Figure) -> go.Figure:
    """Buttons showing assemblies, runs or both. They fade the other kind out (opacity 0, no
    hover) instead of hiding its traces, so a catalogue isolated in the legend stays isolated."""
    kinds = [t.name.split(", ")[-1] for t in fig.data]

    def show(*shown: str) -> dict[str, Any]:
        on = [k in shown or k not in (ASSEMBLIES, RUNS) for k in kinds]
        return {"marker.opacity": [0.65 if o else 0 for o in on], "hoverinfo": ["all" if o else "skip" for o in on]}

    buttons = [
        {"label": "Both", "method": "restyle", "args": [show(ASSEMBLIES, RUNS)]},
        {"label": "Assemblies", "method": "restyle", "args": [show(ASSEMBLIES)]},
        {"label": "Runs", "method": "restyle", "args": [show(RUNS)]},
    ]
    return fig.update_layout(
        updatemenus=[
            {
                "type": "buttons",
                "direction": "right",
                "buttons": buttons,
                "x": 0,
                "y": 0,
                "xanchor": "left",
                "yanchor": "bottom",
            }
        ]
    )


def donuts(rows: pd.DataFrame, records: pd.DataFrame) -> go.Figure:
    """One donut per existing MGnify catalogue: its MAGs now, and the genomes it could gain."""
    existing = rows[rows["catalogue"].notna()].reset_index(drop=True)
    n = max(len(existing), 1)
    fig = make_subplots(
        rows=1,
        cols=n,
        specs=[[{"type": "domain"}] * n],
        subplot_titles=[
            f"{r.biome}<br><sup>{r.catalogue} · updated {r.catalogue_updated}</sup>" for r in existing.itertuples()
        ],
    )
    for i, r in enumerate(existing.itertuples(), 1):
        add = records[
            (records["biome"] == r.biome)
            & records["kind"].isin([SUBMITTED, BY_MGNIFY])
            & (records["in_catalogue"] != "True")
        ]["kind"].value_counts()
        slices = {
            "in catalogue": int(r.catalogue_mags or 0),
            SUBMITTED: int(add.get(SUBMITTED, 0)),
            BY_MGNIFY: int(add.get(BY_MGNIFY, 0)),
        }
        values = list(slices.values())
        shown = {k: v for k, v in slices.items() if v}  # an empty slice would still print its 0
        fig.add_trace(
            go.Pie(
                labels=list(shown),
                values=list(shown.values()),
                hole=0.62,
                sort=False,
                texttemplate="%{value:,}",
                textposition="outside",
                marker_colors=[COLORS[k] for k in shown],
                name=r.biome,
                hovertemplate="%{label}: %{value:,}<extra></extra>",
            ),
            row=1,
            col=i,
        )
        domain = fig.data[-1].domain.x
        fig.add_annotation(
            text=f"<b>+{sum(values[1:]):,}</b><br>MAGs",
            showarrow=False,
            font={"size": 14},
            x=sum(domain) / 2,
            y=0.5,
            xref="paper",
            yref="paper",
        )
    for title in fig.layout.annotations[: len(existing)]:  # clear of the outside slice labels
        title.yshift = 24
    return fig.update_layout(
        FIGURE_LAYOUT
        | {
            "height": 380,
            "margin": {"l": 30, "r": 30, "t": 110, "b": 20},
            "legend": {"orientation": "h", "y": -0.05, "x": 0.5, "xanchor": "center"},
        }
    )


def route_hover(t: pd.DataFrame) -> list[pd.Series]:
    """Hover text of the three routes: the records behind each sample count."""

    def fmt(col: str) -> pd.Series:
        return t[col].map("{:,}".format)

    return [
        fmt("to_add") + " MAGs",
        fmt("to_bin") + " assemblies from " + fmt("samples_bin") + " samples",
        fmt("to_assemble") + " runs from " + fmt("samples_assemble") + " samples",
    ]


def routes_figure(total: pd.Series) -> go.Figure:
    """Where new genomes would come from, for one catalogue."""
    data = pd.DataFrame(
        {
            "route": [f"{k}<br><sup>{v}</sup>" for k, v in ROUTES.items()],
            "n": [total["to_add"], total["samples_bin"], total["samples_assemble"]],
            "detail": [h.iloc[0] for h in route_hover(total.to_frame().T)],
        }
    )
    fig = figure(
        data,
        px.bar,
        x="n",
        y="route",
        orientation="h",
        text_auto=",.0f",
        custom_data=["detail"],
        labels={"n": "", "route": ""},
        layout={"height": 230},
    )
    return fig.update_traces(
        marker_color=[COLORS[SUBMITTED], COLORS[ASSEMBLIES], COLORS[RUNS]],
        hovertemplate="%{customdata[0]}<extra></extra>",
    ).update_yaxes(autorange="reversed")


# --- overview figures (all catalogues) ---


def overview_routes(table_rows: pd.DataFrame) -> go.Figure:
    t = table_rows.set_index("biome")
    columns = [
        "MAGs to add<br><sup>catalogue pipeline</sup>",
        "Samples to bin<br><sup>genomes-generation</sup>",
        "Samples to assemble<br><sup>miassembler</sup>",
    ]
    m = t[["to_add", "samples_bin", "samples_assemble"]].set_axis(columns, axis=1)
    return heatmap(m, "Greens", hover=pd.DataFrame(dict(zip(columns, route_hover(t)))))


def overview_studies(records: pd.DataFrame) -> go.Figure:
    """Projects holding each catalogue's MAGs and the records of each route."""
    genome = records["kind"].isin([SUBMITTED, BY_MGNIFY])
    groups = {
        "MAGs": genome,
        "MAGs to add": genome & (records["in_catalogue"] != "True"),
        "Samples to bin": (records["kind"] == ASSEMBLIES) & (records["route"] == "bin"),
        "Samples to assemble": (records["kind"] == RUNS) & (records["route"] == "assemble"),
    }
    data = pd.DataFrame({k: records[m].groupby("biome")["study_accession"].nunique() for k, m in groups.items()})
    return heatmap(data.reindex(sorted(records["biome"].unique())), "Blues")


# --- catalogue figures (one catalogue) ---


def catalogue_records(taxa: pd.DataFrame) -> go.Figure:
    kinds = ["mags", "samples_bin", "samples_assemble"]
    data = (
        taxa.melt(
            id_vars=["metagenome_taxon"],
            value_vars=kinds,
            var_name="kind",
            value_name="n",
        )
        .replace({"kind": {k: COLUMNS[k][0] for k in COLUMNS}})
        .dropna(subset=["n"])
    )
    fig = figure(
        data,
        px.bar,
        x="n",
        y="metagenome_taxon",
        color="kind",
        facet_col="kind",
        orientation="h",
        text_auto=",.0f",
        color_discrete_map=COLORS,
        labels={"metagenome_taxon": "", "n": ""},
        category_orders={"kind": [COLUMNS[k][0] for k in kinds]},
        layout={"height": 120 + 40 * taxa["metagenome_taxon"].nunique()},
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig.update_xaxes(matches=None, showticklabels=False).update_layout(showlegend=False)


def timeline_figure(records: pd.DataFrame, cutoff: Any) -> go.Figure:
    data = (
        records.dropna(subset=["first_public"])
        .groupby(["kind", pd.Grouper(key="first_public", freq="MS")])
        .size()
        .groupby(level="kind")
        .cumsum()
        .rename("records")
        .reset_index()
    )
    fig = figure(
        data,
        px.line,
        x="first_public",
        y="records",
        color="kind",
        color_discrete_map=COLORS,
        category_orders={"kind": KINDS},
        labels={"first_public": "", "records": "records (cumulative)"},
    )
    if pd.notna(cutoff) and not data.empty:
        fig.add_vline(
            x=pd.Timestamp(cutoff).timestamp() * 1000,
            line_dash="dot",
            annotation_text="catalogue update",
            annotation_position="top left",
        )
    return fig


def top_studies_figure(studies: pd.DataFrame) -> go.Figure:
    """The projects with most records on a route, from studies.tsv."""
    routes = {
        "mags_to_add": "MAGs to add",
        "assemblies_to_bin": "Assemblies to bin",
        "runs_to_assemble": "Runs to assemble",
    }
    top = studies.assign(total=studies[list(routes)].sum(axis=1)).nlargest(TOP_STUDIES, "total")
    label = top["study_accession"] + np.where(top["from_mgnify"] == "True", " (MGnify)", "")
    data = top.assign(study_accession=label).melt(
        id_vars=["study_accession"], value_vars=list(routes), var_name="route", value_name="n"
    )
    data = data[data["n"] > 0].replace({"route": routes})
    return figure(
        data,
        px.bar,
        x="n",
        y="study_accession",
        color="route",
        orientation="h",
        color_discrete_map=COLORS,
        category_orders={"study_accession": list(label), "route": list(routes.values())},
        labels={"n": "records", "study_accession": "", "route": ""},
        layout={"height": 120 + 24 * len(top)},
    )


# --- layout ---


def card(title: str, *objects: Any, tooltip: str = "") -> pn.Card:
    """A titled card; `tooltip` adds an info icon next to the title."""
    header = pn.Row(pn.pane.Markdown(title, margin=(0, 4)), pn.widgets.TooltipIcon(value=tooltip)) if tooltip else None
    return pn.Card(
        *[
            (
                pn.pane.Plotly(o, sizing_mode="stretch_width", config={"displaylogo": False})
                if isinstance(o, go.Figure)
                else o
            )
            for o in objects
        ],
        title=title,
        header=header,
        collapsible=False,
        sizing_mode="stretch_width",
        margin=8,
    )


def grid(*cards: pn.Card) -> pn.GridBox:
    return pn.GridBox(*cards, ncols=2, sizing_mode="stretch_width")


def table(frame: pd.DataFrame, columns: list[str], title: str) -> pn.Card:
    """A compact, width-fitting table: short titles with tooltips, thousands separators."""
    cols = [c for c in columns if c in frame and frame[c].notna().any()]
    data = frame[cols].reset_index(drop=True).astype(object)
    data = data.where(data.notna(), None)
    numeric = [c for c in cols if pd.api.types.is_numeric_dtype(frame[c]) and c != "tax_id"]
    return pn.Card(
        pn.widgets.Tabulator(
            data,
            disabled=True,
            show_index=False,
            layout="fit_columns",
            sizing_mode="stretch_width",
            theme="simple",
            frozen_columns=[cols[0]],
            titles={c: COLUMNS[c][0] for c in cols},
            header_tooltips={c: COLUMNS[c][1] for c in cols},
            widths={cols[0]: 210},
            configuration={"columnDefaults": {"headerWordWrap": True}},
            formatters={c: {"type": "money", "thousand": ",", "precision": 2 if c == "tbases" else 0} for c in numeric},
            text_align={c: "right" for c in numeric},
            header_align="center",
            stylesheets=[TABLE_CSS],
        ),
        title=title,
        collapsible=False,
        sizing_mode="stretch_width",
        margin=8,
    )


def kpis(values: list[tuple[str, Any, str]]) -> pn.FlexBox:
    return pn.FlexBox(
        *[
            pn.indicators.Number(
                label=name,
                value=0 if pd.isna(v) else v,
                format=fmt,
                default_color=ACCENT,
                font_size="24pt",
                title_size="10.5pt",
                margin=(4, 20),
            )
            for name, v, fmt in values
        ],
        margin=(0, 8),
    )


def note(text: str) -> pn.pane.Markdown:
    return pn.pane.Markdown(text, sizing_mode="stretch_width", margin=(0, 12), styles={"font-size": "13px"})


def overview(survey: pd.DataFrame, records: pd.DataFrame, colors: dict[str, str]) -> pn.Column:
    rows = survey[survey["metagenome_taxon"] == "TOTAL"]
    t = totals(rows)
    existing = rows["catalogue"].notna().sum()
    parts = [
        note(ABOUT),
        kpis(
            [
                ("Catalogues", len(rows), "{value}"),
                ("Updates of MGnify catalogues", existing, "{value}"),
                ("New catalogues", len(rows) - existing, "{value}"),
                ("Prokaryote", (rows["kingdom"] == "prokaryotes").sum(), "{value}"),
                ("Eukaryote", (rows["kingdom"] == "eukaryotes").sum(), "{value}"),
            ]
        ),
        table(
            t.assign(biome=t["biome"] + np.where(rows["catalogue"].notna().to_numpy(), " *", "")),
            [
                "biome",
                "kingdom",
                "mags",
                "to_add",
                "samples_bin",
                "samples_assemble",
                "tbases",
            ],
            "Totals per catalogue (* already in MGnify: an update)",
        ),
    ]
    if existing:
        parts.append(
            card(
                "Existing MGnify catalogues: MAGs now, and MAGs not yet in them",
                note(
                    "Grey: MAGs in the current catalogue. Blue: MAGs found by the survey that the "
                    "catalogue doesn't have yet (from submitters, and produced by MGnify)."
                ),
                donuts(rows, records),
            )
        )
    parts += [
        grid(
            card("Where new genomes would come from", overview_routes(t)),
            card("Projects per catalogue", overview_studies(records), tooltip=PROJECTS),
        ),
        grid(card("Reads per run (box: 5th–95th percentile)", depth_figure(records, "biome"))),
        card(
            "Assembly and run locations (1° grid)",
            note(
                "Circles: assemblies; diamonds: runs ≥1M (MAGs have no coordinates in ENA). The buttons "
                "show one kind or both. Double-click a catalogue in the legend to show only it; double-click again to show all."
            ),
            map_figure(records, "biome", colors),
        ),
    ]
    return pn.Column(*parts, sizing_mode="stretch_width")


def headline(total: pd.Series) -> str:
    """Markdown header of one catalogue view, from its TOTAL row."""
    head = f"## {total['biome']}\n**{total['kingdom']}** · "
    if pd.isna(total["catalogue"]):
        return head + "new catalogue (MGnify has none for this biome and kingdom)"
    return (
        head + f"update of MGnify **{total['catalogue']}**: {total['catalogue_species']:,} species, "
        f"{total['catalogue_genomes']:,} genomes ({total['catalogue_mags']:,} MAGs), last updated "
        f"{total['catalogue_updated']}"
    )


def catalogue(rows: pd.DataFrame, records: pd.DataFrame, studies: pd.DataFrame) -> pn.Column:
    raw = rows[rows["metagenome_taxon"] == "TOTAL"].iloc[0]
    t = totals(rows)
    total = t[t["metagenome_taxon"] == "TOTAL"].iloc[0]
    parts = [
        pn.pane.Markdown(headline(raw), sizing_mode="stretch_width", margin=(0, 16)),
        kpis(
            [
                ("MAGs", total["mags"], "{value:,}"),
                ("MAGs to add", total["to_add"], "{value:,}"),
                ("Samples to bin", total["samples_bin"], "{value:,}"),
                ("Samples to assemble", total["samples_assemble"], "{value:,}"),
                ("Runs Tb", total["tbases"], "{value:,.2f}"),
            ]
        ),
    ]
    top = [card("Where new genomes would come from", routes_figure(total))]
    if pd.notna(raw["catalogue"]):
        top.insert(
            0,
            card(
                "MAGs now, and MAGs not yet in the catalogue",
                donuts(rows[rows["metagenome_taxon"] == "TOTAL"], records),
            ),
        )
    parts += [
        grid(*top),
        table(
            t,
            [
                "metagenome_taxon",
                "tax_id",
                "mags",
                "to_add",
                "samples_bin",
                "samples_assemble",
                "tbases",
            ],
            "Counts per metagenome taxon",
        ),
        grid(
            card("MAGs and samples per metagenome taxon", catalogue_records(t[t["metagenome_taxon"] != "TOTAL"])),
            card("Cumulative records by first-public month", timeline_figure(records, raw["catalogue_updated"])),
            card(f"Largest {TOP_STUDIES} projects", top_studies_figure(studies), tooltip=PROJECTS),
            card("Reads per run (box: 5th–95th percentile)", depth_figure(records, "metagenome_taxon")),
        ),
        card(
            "Assembly and run locations (1° grid)",
            note(
                "Circles: assemblies; diamonds: runs ≥1M (MAGs have no coordinates in ENA). "
                "The buttons show one kind or both."
            ),
            map_figure(records, "kind", COLORS),
        ),
    ]
    return pn.Column(*parts, sizing_mode="stretch_width")


def write(out: str | Path) -> None:
    """Write the dashboard `report.html` into the survey folder `out`. Panel templates can't be
    embedded, so it is a plain layout saved with every catalogue view precomputed."""
    out = Path(out)
    pn.extension("plotly", "tabulator", design="bootstrap")
    survey, records = load(out)
    studies = projects(out)
    biomes = sorted(survey["biome"].unique())
    cycle = px.colors.qualitative.Safe + px.colors.qualitative.Vivid
    colors = {b: cycle[i % len(cycle)] for i, b in enumerate(biomes)}
    pick = pn.widgets.Select(label="Catalogue", options=biomes, sizing_mode="stretch_width")

    def view(biome: str) -> pn.Column:
        return catalogue(
            survey[survey["biome"] == biome], records[records["biome"] == biome], studies[studies["biome"] == biome]
        )

    header = pn.pane.Markdown(
        "# biome-finder survey",
        sizing_mode="stretch_width",
        margin=(0, 16),
        styles={"background": ACCENT, "color": "white", "padding": "0 24px", "border-radius": "0 0 6px 6px"},
    )
    page = pn.Column(
        header,
        pn.Tabs(
            ("All catalogues", overview(survey, records, colors)),
            ("One catalogue", pn.Column(pn.Row(pick, max_width=420, margin=(8, 16)), pn.bind(view, pick))),
            dynamic=False,
            sizing_mode="stretch_width",
            margin=(0, 16),
        ),
        sizing_mode="stretch_width",
        max_width=1440,
        styles={"margin": "0 auto"},
    )
    page.save(out / "report.html", title="biome-finder survey", resources=INLINE, embed=True, max_opts=len(biomes))

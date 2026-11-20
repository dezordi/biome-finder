"""Survey dashboard: one self-contained `report.html`, built with Panel and plotly.

Drawn only from the tables in a survey folder, so it can be redrawn without querying ENA.
Two tabs: an overview of all catalogues, and one catalogue picked from a dropdown. Every
catalogue view is precomputed into the page (`save(embed=True)`), so the file works offline.

Everything shown is for the catalogue's kingdom: its genomes, and only the assemblies and runs
with evidence for that kingdom (`kingdom_match` in records.tsv). The overview compares many
catalogues, so it uses heatmaps with the count in each cell; a catalogue view uses bars.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import panel as pn
import plotly.express as px
import plotly.graph_objects as go
from bokeh.resources import INLINE
from plotly.subplots import make_subplots

SUBMITTED, BY_MGNIFY = "MAGs (submitters)", "MAGs (MGnify)"
ASSEMBLIES, RUNS = "Assemblies", "Runs ≥1M"
KINDS = [SUBMITTED, BY_MGNIFY, ASSEMBLIES, RUNS]
ROUTES = {"Genomes to add": "→ genomes-catalogue-pipeline", "Assemblies to bin": "→ genomes-generation",
          "Runs to assemble": "→ miassembler, then genomes-generation"}
LONG_READS = {"OXFORD_NANOPORE", "PACBIO_SMRT"}
TOP_STUDIES = 15
ACCENT = "#0072B2"
COLORS = {SUBMITTED: "#0072B2", BY_MGNIFY: "#56B4E9", ASSEMBLIES: "#E69F00", RUNS: "#009E73",
          "short reads": "#009E73", "long reads": "#CC79A7", "in catalogue": "#BBBBBB",
          "Bacteria": "#0072B2", "Archaea": "#D55E00", "Eukaryota": "#009E73"}
# column: (title, tooltip)
COLUMNS = {
    "biome": ("Catalogue", "The biome config: one catalogue per biome and kingdom, as in MGnify"),
    "metagenome_taxon": ("Metagenome taxon", "ENA metagenome taxon the records were found under"),
    "tax_id": ("Tax ID", "NCBI tax ID of the metagenome taxon"),
    "kingdom": ("Kingdom", "prokaryotes (Bacteria + Archaea) or eukaryotes"),
    "samples": ("Samples", "Distinct biological samples behind the kingdom's genomes (their source samples), "
                           "assemblies and runs"),
    "submitted": (SUBMITTED, "MAGs of the catalogue's kingdom deposited in ENA by their authors"),
    "by_mgnify": (BY_MGNIFY, "Genomes MGnify produced from this biome's public data. ENA stores most of them as "
                             "'binned metagenome' records, a few as MAGs"),
    "assemblies": (ASSEMBLIES, "Metagenome assemblies with evidence for the kingdom: a genome of that kingdom from the "
                               "same sample, or a project that names it (prokaryotes: all but eukaryote-only projects)"),
    "runs": (RUNS, "Metagenomic WGS runs with at least 1 million reads and evidence for the kingdom"),
    "gbases": ("Gbases", "Gigabases in those runs"),
    "new_runs": ("New runs", "Of those runs, public after the MGnify catalogue's last update (blank: new catalogue)"),
    "to_add": ("Genomes to add", "MAGs (both kinds) not in the existing MGnify catalogue (all of them, for a new "
                                 "catalogue): ready for genomes-catalogue-pipeline"),
    "to_bin": ("Assemblies to bin", "Assemblies whose sample has no genome yet: input for genomes-generation"),
    "to_assemble": ("Runs to assemble", "Runs ≥1M whose sample has neither a genome nor an assembly: input for "
                                        "miassembler, then genomes-generation"),
}
FIGURE_LAYOUT = {"template": "plotly_white", "height": 380, "margin": {"l": 10, "r": 10, "t": 36, "b": 10},
                 "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0, "title": None},
                 "font": {"size": 12}}
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
*genomes to add* are ready for the catalogue; *assemblies to bin* need binning first; *runs to
assemble* need assembly and binning. Many genomes to add means the catalogue can grow now; many
assemblies or runs means it would pay to generate new MAGs.
"""


def load(out):
    """survey.tsv plus every catalogue's records.tsv, with a `biome` column and a `kind`
    (MAGs from submitters or MGnify, assemblies, runs ≥1M) for every record counted for the
    catalogue's kingdom; everything else is dropped."""
    survey = pd.read_csv(out / "survey.tsv", sep="\t", dtype={"catalogue_updated": str, "tax_id": str},
                         dtype_backend="numpy_nullable")
    records = pd.concat([pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False).assign(biome=p.parent.name)
                         for p in sorted((out / "biomes").glob("*/records.tsv"))], ignore_index=True)
    records["first_public"] = pd.to_datetime(records["first_public"], errors="coerce")
    for col in ["read_count", "base_count", "lat", "lon"]:
        records[col] = pd.to_numeric(records[col], errors="coerce")
    genome = records["data_type"].isin(["mag", "bin"])
    by_mgnify = (records["data_type"] == "bin") | (records["from_mgnify"] == "True")
    matched = records["kingdom_match"] == "True"
    records["kind"] = np.select(
        [genome & by_mgnify, genome, (records["data_type"] == "contigs") & matched,
         (records["data_type"] == "reads") & matched & (records["read_count"] >= 1e6)],
        [BY_MGNIFY, SUBMITTED, ASSEMBLIES, RUNS], default="")
    return survey, records[records["kind"] != ""].copy()


def samples(records):
    """Distinct samples behind the kingdom records, per (biome, metagenome taxon) plus TOTAL.
    Genomes count their source samples (`derived from`), so a co-assembly counts all of them."""
    src = records["source_samples"].where(records["kind"].isin([SUBMITTED, BY_MGNIFY]) & (records["source_samples"] != ""),
                                          records["sample_accession"])
    # ponytail: primary and secondary IDs (SAMEA vs ERS) of one sample can both appear; map them if it matters
    s = records.assign(sample=src.str.split(";")).explode("sample")
    s = s[s["sample"] != ""]
    per_taxon = s.groupby(["biome", "metagenome_taxon"])["sample"].nunique()
    total = s.groupby("biome")["sample"].nunique()
    total.index = pd.MultiIndex.from_product([total.index, ["TOTAL"]])
    return pd.concat([per_taxon, total])


def totals(rows, records):
    """survey.tsv rows as the dashboard's kingdom-only columns."""
    n = rows.fillna({c: 0 for c in rows.columns if pd.api.types.is_numeric_dtype(rows[c])})
    return pd.DataFrame({
        "biome": rows["biome"], "metagenome_taxon": rows["metagenome_taxon"], "tax_id": rows["tax_id"],
        "kingdom": rows["kingdom"],
        "samples": [samples(records).get((b, x), 0) for b, x in zip(rows["biome"], rows["metagenome_taxon"])],
        "submitted": n["mags"] - n["mags_mgnify"], "by_mgnify": n["mags_mgnify"] + n["mgnify_bins"],
        "assemblies": rows["contigs_kingdom"], "runs": rows["runs_1m_kingdom"], "gbases": rows["gbases_1m_kingdom"],
        "new_runs": rows["runs_1m_kingdom_new"], "to_add": n["mags_to_add"] + n["mgnify_bins_to_add"],
        "to_bin": rows["contigs_to_bin"], "to_assemble": rows["runs_to_assemble"],
    }).reset_index(drop=True)


# --- figure helpers ---

def figure(frame, plot, layout=None, **kw):
    """A styled plotly express figure, or a 'no data' placeholder for an empty frame."""
    if frame.empty:
        fig = go.Figure()
        fig.add_annotation(text="no data", showarrow=False, font={"size": 16, "color": "grey"})
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
    else:
        fig = plot(frame, **kw)
    return fig.update_layout(FIGURE_LAYOUT | (layout or {}))


def heatmap(matrix, colorscale="Blues"):
    """Counts as a heatmap: log-scaled color, the count written in each cell, height to fit."""
    matrix = matrix.astype(float)
    text = matrix.map(lambda v: "" if pd.isna(v) or v == 0 else f"{v:,.0f}")
    fig = go.Figure(go.Heatmap(
        z=np.log10(matrix.where(matrix > 0)), x=[str(c) for c in matrix.columns], y=list(matrix.index),
        text=text.values, texttemplate="%{text}", colorscale=colorscale, showscale=False, xgap=2, ygap=2,
        hovertemplate="%{y} · %{x}: %{text}<extra></extra>"))
    fig.update_yaxes(autorange="reversed", ticksuffix="  ")
    fig.update_xaxes(side="top", type="category", tickangle=0)
    return fig.update_layout(FIGURE_LAYOUT | {"height": 70 + 26 * len(matrix),
                                              "margin": {"l": 10, "r": 10, "t": 50, "b": 10}})


def depth_figure(records, by):
    """Reads per run as boxes from precomputed percentiles (5/25/50/75/95)."""
    runs = records[records["kind"] == RUNS]
    runs = runs.assign(platform=runs["instrument_platform"].isin(LONG_READS).map(
        {True: "long reads", False: "short reads"}))
    q = runs.groupby(["platform", by])["read_count"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
    fig = go.Figure([go.Box(name=platform, y=list(g.index.get_level_values(by)), lowerfence=g[0.05],
                            q1=g[0.25], median=g[0.5], q3=g[0.75], upperfence=g[0.95], orientation="h",
                            marker_color=COLORS[platform])
                     for platform, g in q.groupby(level="platform")])
    n = max(runs[by].nunique(), 1)
    return (fig.update_layout(FIGURE_LAYOUT | {"height": max(260, 80 + 34 * n)}, boxmode="group")
            .update_xaxes(type="log", title="reads per run").update_yaxes(categoryorder="category descending"))


def map_figure(records, color, colors):
    data = (records[records["kind"].isin([ASSEMBLIES, RUNS])].dropna(subset=["lat", "lon"])
            .assign(lat=lambda d: d["lat"].round(), lon=lambda d: d["lon"].round())
            .groupby([color, "lat", "lon"], as_index=False)
            .agg(records=("accession", "size"), studies=("study_accession", "nunique")))
    fig = figure(data, px.scatter_geo, layout={"height": 520}, lat="lat", lon="lon", size="records",
                 color=color, color_discrete_map=colors, hover_data=["records", "studies"], size_max=22,
                 opacity=0.65, projection="natural earth")
    return fig.update_geos(showcountries=True, countrycolor="#d0d0d0", showland=True, landcolor="#f4f4f4",
                           showframe=False)


def donuts(rows, records):
    """One donut per existing MGnify catalogue: its MAGs now, and the genomes it could gain."""
    existing = rows[rows["catalogue"].notna()].reset_index(drop=True)
    n = max(len(existing), 1)
    fig = make_subplots(rows=1, cols=n, specs=[[{"type": "domain"}] * n],
                        subplot_titles=[f"{r.biome}<br><sup>{r.catalogue} · updated {r.catalogue_updated}</sup>"
                                        for r in existing.itertuples()])
    for i, r in enumerate(existing.itertuples(), 1):
        add = records[(records["biome"] == r.biome) & records["kind"].isin([SUBMITTED, BY_MGNIFY])
                      & (records["in_catalogue"] != "True")]["kind"].value_counts()
        slices = {"in catalogue": int(r.catalogue_mags or 0), SUBMITTED: int(add.get(SUBMITTED, 0)),
                  BY_MGNIFY: int(add.get(BY_MGNIFY, 0))}
        values = list(slices.values())
        shown = {k: v for k, v in slices.items() if v}  # an empty slice would still print its 0
        fig.add_trace(go.Pie(labels=list(shown), values=list(shown.values()), hole=0.62, sort=False,
                             texttemplate="%{value:,}", textposition="outside",
                             marker_colors=[COLORS[k] for k in shown], name=r.biome,
                             hovertemplate="%{label}: %{value:,}<extra></extra>"), row=1, col=i)
        domain = fig.data[-1].domain.x
        fig.add_annotation(text=f"<b>+{sum(values[1:]):,}</b><br>genomes", showarrow=False, font={"size": 14},
                           x=sum(domain) / 2, y=0.5, xref="paper", yref="paper")
    for title in fig.layout.annotations[:len(existing)]:  # clear of the outside slice labels
        title.yshift = 24
    return fig.update_layout(FIGURE_LAYOUT | {"height": 380, "margin": {"l": 30, "r": 30, "t": 110, "b": 20},
                                              "legend": {"orientation": "h", "y": -0.05, "x": 0.5,
                                                         "xanchor": "center"}})


def routes_figure(total):
    """Where new genomes would come from, for one catalogue."""
    data = pd.DataFrame({"route": [f"{k}<br><sup>{v}</sup>" for k, v in ROUTES.items()],
                         "n": [total["to_add"], total["to_bin"], total["to_assemble"]]})
    fig = figure(data, px.bar, x="n", y="route", orientation="h", text_auto=",.0f",
                 labels={"n": "", "route": ""}, layout={"height": 230})
    return fig.update_traces(marker_color=[COLORS[SUBMITTED], COLORS[ASSEMBLIES], COLORS[RUNS]]).update_yaxes(
        autorange="reversed")


# --- overview figures (all catalogues) ---

def overview_routes(table_rows):
    m = table_rows.set_index("biome")[["to_add", "to_bin", "to_assemble"]]
    return heatmap(m.rename(columns={"to_add": "Genomes to add<br><sup>catalogue pipeline</sup>",
                                     "to_bin": "Assemblies to bin<br><sup>genomes-generation</sup>",
                                     "to_assemble": "Runs to assemble<br><sup>miassembler</sup>"}), "Greens")


def overview_studies(records):
    data = records.groupby(["biome", "kind"])["study_accession"].nunique().unstack()
    return heatmap(data.reindex(columns=[k for k in KINDS if k in data]), "Blues")


def overview_years(records):
    return heatmap(records.dropna(subset=["first_public"]).groupby(
        ["biome", records["first_public"].dt.year]).size().unstack().sort_index(axis=1), "Purples")


def overview_concentration(records):
    per_study = records.groupby(["biome", "study_accession"]).size().rename("n").reset_index()
    top = (per_study.sort_values("n", ascending=False).groupby("biome")
           .agg(top=("n", lambda r: r.head(3).sum()), total=("n", "sum")))
    data = (top["top"] / top["total"]).rename("share").reset_index().sort_values("share")
    fig = figure(data, px.bar, x="share", y="biome", orientation="h", text_auto=".0%",
                 labels={"share": "share of records in the 3 largest studies", "biome": ""},
                 layout={"height": 70 + 26 * len(data)})
    return fig.update_traces(marker_color=ACCENT).update_xaxes(tickformat=".0%", range=[0, 1])


def overview_groups(records):
    genomes = records[records["kind"].isin([SUBMITTED, BY_MGNIFY])]
    if genomes.empty:
        return figure(genomes, px.bar)
    return heatmap(genomes.groupby(["biome", "scope_taxon"]).size().unstack(), "Oranges")


# --- catalogue figures (one catalogue) ---

def catalogue_records(taxa):
    data = (taxa.melt(id_vars=["metagenome_taxon"], value_vars=["submitted", "by_mgnify", "assemblies", "runs"],
                      var_name="kind", value_name="n")
            .replace({"kind": {k: COLUMNS[k][0] for k in COLUMNS}}).dropna(subset=["n"]))
    fig = figure(data, px.bar, x="n", y="metagenome_taxon", color="kind", facet_col="kind", orientation="h",
                 text_auto=",.0f", color_discrete_map=COLORS, labels={"metagenome_taxon": "", "n": ""},
                 category_orders={"kind": KINDS}, layout={"height": 120 + 40 * taxa["metagenome_taxon"].nunique()})
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig.update_xaxes(matches=None, showticklabels=False).update_layout(showlegend=False)


def timeline_figure(records, cutoff):
    data = (records.dropna(subset=["first_public"])
            .groupby(["kind", pd.Grouper(key="first_public", freq="MS")]).size()
            .groupby(level="kind").cumsum().rename("records").reset_index())
    fig = figure(data, px.line, x="first_public", y="records", color="kind", color_discrete_map=COLORS,
                 category_orders={"kind": KINDS}, labels={"first_public": "", "records": "records (cumulative)"})
    if pd.notna(cutoff) and not data.empty:
        fig.add_vline(x=pd.Timestamp(cutoff).timestamp() * 1000, line_dash="dot",
                      annotation_text="catalogue update", annotation_position="top left")
    return fig


def top_studies_figure(records):
    per = records.groupby(["study_accession", "kind"]).size().rename("n").reset_index()
    mgnify = set(records.loc[records["from_mgnify"] == "True", "study_accession"])
    order = per.groupby("study_accession")["n"].sum().nlargest(TOP_STUDIES).index
    data = per[per["study_accession"].isin(order)]
    label = {s: f"{s} (MGnify)" if s in mgnify else s for s in order}
    data = data.assign(study_accession=data["study_accession"].map(label))
    return figure(data, px.bar, x="n", y="study_accession", color="kind", orientation="h",
                  color_discrete_map=COLORS, category_orders={"study_accession": [label[s] for s in order],
                                                              "kind": KINDS},
                  labels={"n": "records", "study_accession": "", "kind": ""}, layout={"height": 120 + 24 * len(order)})


def concentration_figure(records):
    per = records.groupby("study_accession").size().sort_values(ascending=False)
    data = pd.DataFrame({"rank": range(1, len(per) + 1), "share": per.cumsum() / per.sum(),
                         "study": per.index, "records": per.values})
    fig = figure(data, px.line, x="rank", y="share", log_x=True, markers=len(data) < 30,
                 hover_data=["study", "records"],
                 labels={"rank": "studies, largest first", "share": "share of records"})
    return fig.update_traces(line_color=ACCENT).update_yaxes(tickformat=".0%", range=[0, 1.02])


def catalogue_groups(records):
    data = (records[records["kind"].isin([SUBMITTED, BY_MGNIFY])]
            .groupby(["scope_taxon", "kind"], as_index=False).size())
    return figure(data, px.bar, x="size", y="scope_taxon", color="kind", orientation="h", text_auto=",.0f",
                  color_discrete_map=COLORS, category_orders={"kind": KINDS},
                  labels={"size": "genomes", "scope_taxon": "", "kind": ""}, layout={"height": 220})


# --- layout ---

def card(title, *objects):
    return pn.Card(*[pn.pane.Plotly(o, sizing_mode="stretch_width", config={"displaylogo": False})
                     if isinstance(o, go.Figure) else o for o in objects],
                   title=title, collapsible=False, sizing_mode="stretch_width", margin=8)


def grid(*cards):
    return pn.GridBox(*cards, ncols=2, sizing_mode="stretch_width")


def table(frame, columns, title):
    """A compact, width-fitting table: short titles with tooltips, thousands separators."""
    cols = [c for c in columns if c in frame and frame[c].notna().any()]
    data = frame[cols].reset_index(drop=True).astype(object)
    data = data.where(data.notna(), None)
    numeric = [c for c in cols if pd.api.types.is_numeric_dtype(frame[c]) and c != "tax_id"]
    return pn.Card(pn.widgets.Tabulator(
        data, disabled=True, show_index=False, layout="fit_columns", sizing_mode="stretch_width",
        theme="simple", frozen_columns=[cols[0]], titles={c: COLUMNS[c][0] for c in cols},
        header_tooltips={c: COLUMNS[c][1] for c in cols}, widths={cols[0]: 210},
        configuration={"columnDefaults": {"headerWordWrap": True}},
        formatters={c: {"type": "money", "thousand": ",", "precision": 1 if c == "gbases" else 0} for c in numeric},
        text_align={c: "right" for c in numeric}, header_align="center", stylesheets=[TABLE_CSS]),
        title=title, collapsible=False, sizing_mode="stretch_width", margin=8)


def kpis(values):
    return pn.FlexBox(*[pn.indicators.Number(label=name, value=0 if pd.isna(v) else v, format=fmt,
                                             default_color=ACCENT, font_size="24pt", title_size="10.5pt",
                                             margin=(4, 20))
                        for name, v, fmt in values], margin=(0, 8))


def note(text):
    return pn.pane.Markdown(text, sizing_mode="stretch_width", margin=(0, 12), styles={"font-size": "13px"})


def overview(survey, records, colors):
    rows = survey[survey["metagenome_taxon"] == "TOTAL"]
    t = totals(rows, records)
    existing = rows["catalogue"].notna().sum()
    parts = [
        note(ABOUT),
        kpis([("Catalogues", len(rows), "{value}"), ("Updates of MGnify catalogues", existing, "{value}"),
              ("New catalogues", len(rows) - existing, "{value}"),
              ("Prokaryote", (rows["kingdom"] == "prokaryotes").sum(), "{value}"),
              ("Eukaryote", (rows["kingdom"] == "eukaryotes").sum(), "{value}")]),
        table(t, ["biome", "kingdom", "samples", "submitted", "by_mgnify", "assemblies", "runs", "gbases",
                  "new_runs", "to_add", "to_bin", "to_assemble"], "Totals per catalogue"),
    ]
    if existing:
        parts.append(card("Existing MGnify catalogues: MAGs now, and genomes not yet in them",
                          note("Grey: MAGs in the current catalogue. Blue: MAGs found by the survey that the "
                               "catalogue doesn't have yet (from submitters, and produced by MGnify)."),
                          donuts(rows, records)))
    parts += [
        grid(card("Where new genomes would come from", overview_routes(t)),
             card("Studies per catalogue", overview_studies(records)),
             card("Genomes per organism group", overview_groups(records)),
             card("Share of records in the 3 largest studies", overview_concentration(records))),
        card("Records first public, per year", overview_years(records)),
        grid(card("Reads per run (box: 5th–95th percentile)", depth_figure(records, "biome"))),
        card("Assembly and run locations (1° grid)", map_figure(records, "biome", colors)),
    ]
    return pn.Column(*parts, sizing_mode="stretch_width")


def headline(total):
    """Markdown header of one catalogue view, from its TOTAL row."""
    head = f"## {total['biome']}\n**{total['kingdom']}** · "
    if pd.isna(total["catalogue"]):
        return head + "new catalogue (MGnify has none for this biome and kingdom)"
    return (head + f"update of MGnify **{total['catalogue']}**: {total['catalogue_species']:,} species, "
            f"{total['catalogue_genomes']:,} genomes ({total['catalogue_mags']:,} MAGs), last updated "
            f"{total['catalogue_updated']}")


def catalogue(rows, records):
    raw = rows[rows["metagenome_taxon"] == "TOTAL"].iloc[0]
    t = totals(rows, records)
    total = t[t["metagenome_taxon"] == "TOTAL"].iloc[0]
    parts = [
        pn.pane.Markdown(headline(raw), sizing_mode="stretch_width", margin=(0, 16)),
        kpis([("Samples", total["samples"], "{value:,}"),
              (SUBMITTED, total["submitted"], "{value:,}"), (BY_MGNIFY, total["by_mgnify"], "{value:,}"),
              (ASSEMBLIES, total["assemblies"], "{value:,}"), (RUNS, total["runs"], "{value:,}"),
              ("Gbases", total["gbases"], "{value:,.0f}"), ("New runs", total["new_runs"], "{value:,}")]),
    ]
    top = [card("Where new genomes would come from", routes_figure(total))]
    if pd.notna(raw["catalogue"]):
        top.insert(0, card("MAGs now, and genomes not yet in the catalogue",
                           donuts(rows[rows["metagenome_taxon"] == "TOTAL"], records)))
    parts += [
        grid(*top),
        table(t, ["metagenome_taxon", "tax_id", "submitted", "by_mgnify", "assemblies", "runs", "gbases", "new_runs",
                  "to_add", "to_bin", "to_assemble"], "Counts per metagenome taxon"),
        grid(card("Records per metagenome taxon", catalogue_records(t[t["metagenome_taxon"] != "TOTAL"])),
             card("Cumulative records by first-public month", timeline_figure(records, raw["catalogue_updated"])),
             card(f"Largest {TOP_STUDIES} studies", top_studies_figure(records)),
             card("Share of records in the largest studies", concentration_figure(records)),
             card("Genomes per organism group", catalogue_groups(records)),
             card("Reads per run (box: 5th–95th percentile)", depth_figure(records, "metagenome_taxon"))),
        card("Assembly and run locations (1° grid)", map_figure(records, "kind", COLORS)),
    ]
    return pn.Column(*parts, sizing_mode="stretch_width")


def write(out):
    """Write the dashboard `report.html` into the survey folder `out`."""
    out = Path(out)
    pn.extension("plotly", "tabulator", design="bootstrap")
    survey, records = load(out)
    biomes = sorted(survey["biome"].unique())
    cycle = px.colors.qualitative.Safe + px.colors.qualitative.Vivid
    colors = {b: cycle[i % len(cycle)] for i, b in enumerate(biomes)}
    pick = pn.widgets.Select(label="Catalogue", options=biomes, sizing_mode="stretch_width")

    def view(biome):
        return catalogue(survey[survey["biome"] == biome], records[records["biome"] == biome])

    header = pn.pane.Markdown("# biome-finder survey", sizing_mode="stretch_width", margin=(0, 16),
                              styles={"background": ACCENT, "color": "white", "padding": "0 24px",
                                      "border-radius": "0 0 6px 6px"})
    page = pn.Column(
        header,
        pn.Tabs(("All catalogues", overview(survey, records, colors)),
                ("One catalogue", pn.Column(pn.Row(pick, max_width=420, margin=(8, 16)), pn.bind(view, pick))),
                dynamic=False, sizing_mode="stretch_width", margin=(0, 16)),
        sizing_mode="stretch_width", max_width=1440, styles={"margin": "0 auto"},
    )
    # Panel templates can't be embedded; a plain layout can, with every catalogue precomputed
    page.save(out / "report.html", title="biome-finder survey", resources=INLINE, embed=True,
              max_opts=len(biomes))

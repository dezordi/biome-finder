"""Writes docs/img/survey-flow.svg: the marine-water-euk config (left) and the survey steps (right).

The numbers come from a survey run; edit STEPS when they change, then:
    python docs/img/survey_flow.py docs/img/survey-flow.svg
"""

import sys
from html import escape

STEPS = [  # (title, detail lines, example lines)
    (
        "Catalogue context",
        ["MGnify API v2 · type must match kingdom · last update = cutoff"],
        ["eukaryotes ✓ · 74 species · “new” = public since 2026-10-27"],
    ),
    (
        "Metagenome taxa → tax IDs",
        ["each name → its ENA tax ID → three searches and one output row"],
        ["marine metagenome → 408172 · seawater metagenome → 1561972", "marine plankton metagenome → 1874687"],
    ),
    (
        "ENA searches, per taxon",
        [
            "MAGs, bins: only the config's kingdom, minus the excluded taxa",
            "assemblies, runs: whole communities; eukaryote subset in step 6",
        ],
        [
            "61 eukaryote MAGs (of 213,712 of any kingdom) + 235 MGnify bins (step 5)",
            "4,616 eukaryote assemblies (of 20,349) · 5,346 runs (of 30,542)",
        ],
    ),
    (
        "Host filter (host-associated biomes only)",
        [
            "runs, assemblies: sample host tax ID = host_tax_id,",
            "or the free-text host contains one of host_keywords",
            "MAGs: derived from → source samples; kept only if all pass",
        ],
        ["not applied: marine water has no host"],
    ),
    (
        "MGnify's own data",
        [
            "records in MGnify's ENA studies (center contains EMG), plus the",
            "bins MGnify made from this biome's samples (binned metagenome)",
        ],
        [
            "38 MAGs + 235 bins, all eukaryote (259 of the catalogue's 294 MAGs)",
            "1,356 eukaryote assemblies (of 6,112 MGnify assemblies)",
        ],
    ),
    (
        "Kingdom evidence (assemblies, runs)",
        [
            "genome: the sample already yielded a genome of the kingdom",
            "text: one of the sample's projects names the kingdom",
            "eukaryotes need one of these; prokaryotes skip only euk-only projects",
        ],
        [
            "eukaryote assemblies: 4,616 (411 genome · 4,205 text) of 20,349",
            "eukaryote runs ≥ 1 M: 3,427 (134 genome · 3,293 text) of 22,741",
        ],
    ),
    (
        "Summarize",
        ["per taxon + TOTAL · runs ≥ 1 M reads · paired / long · new"],
        [
            "eukaryote: 61 MAGs · 235 MGnify bins · 4,616 assemblies",
            "3,427 runs ≥ 1 M (86,629 Gbases) · 12 new since the cutoff",
        ],
    ),
]
YAML = [  # (text, key?, badges)
    ("name:", " marine-water-euk", [7]),
    ("kingdom:", " eukaryotes", [1, 3, 6]),
    ("mgnify_catalogue:", " marine-eukaryotes-v1-0", [1]),
    ("metagenome_taxa:", "", [2]),
    (None, "  - marine metagenome", []),
    (None, "  - seawater metagenome", []),
    (None, "  - marine plankton metagenome", []),
    ("exclude_taxa_extra:", " []", [3]),
    ("host_tax_id:", "", [4]),
    ("host_keywords:", " []", [4]),
]
COMMENTS = [
    "# not used by survey (search, prepare):",
    "# catalogue_biome, mgnify_lineages*,",
    "# envo_*, keywords_*,",
    "# gold_ecosystems, data_types, …",
]
GAP, LINE, X, W = 20, 18, 490, 490


def text(x, y, s, cls, extra=""):
    return f'<text x="{x}" y="{y}" class="{cls}"{extra}>{escape(s).replace("  ", "&#160;&#160;")}</text>'


def badge(cx, cy, n, r=9):
    return (
        f'<circle cx="{cx}" cy="{cy}" r="{r}" class="badge"/>'
        f'<text x="{cx}" y="{cy + 4}" class="badge-text">{n}</text>'
    )


out, y = [], 52
for i, (title, details, examples) in enumerate(STEPS, 1):
    h = 50 + LINE * (len(details) - 1) + 22 + 20 * (len(examples) - 1) + 14
    out.append(f'<rect x="{X}" y="{y}" width="{W}" height="{h}" rx="6" class="box"/>')
    out.append(badge(X + 22, y + 22, i, 11))
    out.append(text(X + 42, y + 27, title, "title"))
    ty = y + 50
    for d in details:
        out.append(text(X + 22, ty, d, "detail"))
        ty += LINE
    ty += 4
    for e in examples:
        out.append(text(X + 22, ty, e, "example"))
        ty += 20
    y += h
    out.append(
        f'<line x1="{X + W // 2}" y1="{y}" x2="{X + W // 2}" y2="{y + GAP - 4}" class="flow" marker-end="url(#arrow)"/>'
    )
    y += GAP
out.append(f'<rect x="{X}" y="{y}" width="{W}" height="62" rx="6" fill="#1f2933"/>')
out.append(
    f'<text x="{X + W // 2}" y="{y + 26}" fill="#ffffff" font-size="13.5px" font-weight="bold" text-anchor="middle">survey/</text>'
)
out.append(
    f'<text x="{X + W // 2}" y="{y + 48}" fill="#dfe6ee" font-size="12.5px" text-anchor="middle">'
    f'{escape("survey.tsv · biomes/<name>/records.tsv, studies.tsv · report.html")}</text>'
)
height = y + 62 + 24

left = []
ly = 82
for key, rest, badges in YAML:
    body = (f'<tspan class="key">{escape(key)}</tspan>' if key else "") + escape(rest).replace("  ", "&#160;&#160;")
    left.append(f'<text x="36" y="{ly}">{body}</text>')
    for k, n in enumerate(reversed(badges)):
        left.append(badge(430 - 22 * k, ly - 5, n))
    ly += 22
ly += 22
for c in COMMENTS:
    left.append(text(36, ly, c, "comment"))
    ly += 22
box_h = ly - 52 - 4

svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="{height}" viewBox="0 0 1000 {height}" font-family="Helvetica, Arial, sans-serif">
  <title>How a biome config drives biome-finder survey (example: marine-water-euk)</title>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#7a8a99"/></marker>
    <style>
      .yaml text:not(.badge-text), text.comment {{ font-family: Menlo, Consolas, "DejaVu Sans Mono", monospace; font-size: 13px; fill: #1f2933; }}
      .key {{ fill: #0b4f8a; font-weight: bold; }}
      text.comment {{ fill: #8a94a0; }}
      .badge {{ fill: #0072B2; }}
      .badge-text {{ fill: #ffffff; font-size: 11px; font-weight: bold; text-anchor: middle; font-family: Helvetica, Arial, sans-serif; }}
      .box {{ fill: #f5f8fb; stroke: #c9d6e3; stroke-width: 1.2; }}
      .title {{ font-size: 15px; font-weight: bold; fill: #1f2933; }}
      .detail {{ font-size: 12.5px; fill: #44515e; }}
      .example {{ font-size: 12.5px; fill: #0072B2; font-style: italic; }}
      .heading {{ font-size: 14px; font-weight: bold; fill: #1f2933; }}
      .flow {{ stroke: #7a8a99; stroke-width: 1.5; }}
    </style>
  </defs>
  <rect width="1000" height="{height}" fill="#ffffff"/>
  <text x="20" y="38" class="heading">biomes/marine-water-euk.yaml</text>
  <rect x="20" y="52" width="430" height="{box_h}" rx="6" fill="#fbfbf8" stroke="#d8d8cf"/>
  <g class="yaml">
    {chr(10).join("    " + l for l in left).strip()}
  </g>
  <text x="490" y="38" class="heading">survey steps</text>
  {chr(10).join("  " + l for l in out).strip()}
</svg>
"""
open(sys.argv[1], "w").write(svg)

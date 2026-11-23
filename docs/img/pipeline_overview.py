"""Writes docs/img/pipeline-overview.svg: the biome-finder steps (left) and the data each one reads (right).

python docs/img/pipeline_overview.py docs/img/pipeline-overview.svg
"""

import sys
from html import escape

# (command, status, what it does, output, data it reads, data it does NOT read)
STAGES = [
    (
        "survey",
        "implemented",
        [
            "Counts the public data per candidate biome, before any curation:",
            "MAGs, MGnify bins, assemblies and runs of the config's kingdom.",
            "Decides which catalogues to build, or which MGnify ones to update.",
        ],
        "survey/  survey.tsv · biomes/<name>/ · report.html",
        [
            "ENA Portal API: MAGs, bins, assemblies, runs, samples, studies",
            "ENA Browser XML: derived from (genome → source samples)",
            "ENA taxonomy: metagenome taxa → tax IDs",
            "MGnify API v1/v2 + FTP: existing catalogues (context only)",
        ],
        ["no papers, no LLM, no NCBI, no GOLD, no composition", "only metadata already public in ENA"],
    ),
    (
        "search",
        "planned",
        [
            "One row per record: data type, resolved taxonomy, biome evidence,",
            "score and tier. Rule-based actions first; papers, LLM and",
            "composition only for records that are still uncertain.",
        ],
        "samples.tsv",
        [
            "ENA Portal / Browser / xref · NCBI E-utilities, FTP, taxdump",
            "MGnify, GOLD, Metalog, mag-provenance",
            "papers: Europe PMC open access + curator PDFs (--papers-dir)",
            "LLM on the retrieved text (--llm), quotes verified",
            "composition (--composition): Trapiche on MGnify taxonomy files",
            "  + our text evidence; IndicPiper on Sandpiper GTDB profiles",
        ],
        [],
    ),
    (
        "prepare",
        "planned",
        [
            "Keeps the curated rows and routes each sample by its most",
            "processed data. Writes status and paths back to samples.tsv.",
        ],
        "catalogue/ · binning/ · assembly/",
        [
            "ENA / NCBI genome FASTAs (vendored fetch_ena.py, fetch_ncbi.py)",
            "ENA FTP URLs + md5 for contigs and reads (not downloaded)",
        ],
        [],
    ),
    (
        "link-assemblies",
        "planned",
        ["After miassembler: builds the genomes-generation inputs", "from its assemblies."],
        "binning/",
        ["miassembler output folder (local)"],
        [],
    ),
]
HUMAN = [  # between stage i and i+1
    "curator: picks 1–2 biomes and refines their YAML (define-biome)",
    "curator: decides tier C / needs_curation rows, adds PDFs, re-runs search",
    "run miassembler on assembly/samplesheet.csv",
]
DOWNSTREAM = "genomes-catalogue-pipeline · genomes-generation · miassembler (outside biome-finder)"
LX, LW, RX, RW, LINE, GAP = 20, 520, 570, 510, 18, 64


def text(x, y, s, cls, extra=""):
    return f'<text x="{x}" y="{y}" class="{cls}"{extra}>{escape(s)}</text>'


out, y = [], 60
for i, (cmd, status, does, output, reads, skips) in enumerate(STAGES):
    right = len(reads) + (len(skips) + 1 if skips else 0)
    h = max(52 + LINE * len(does) + 26, 26 + LINE * right + 6)
    out.append(f'<rect x="{LX}" y="{y}" width="{LW}" height="{h}" rx="6" class="box"/>')
    out.append(
        f'<circle cx="{LX + 22}" cy="{y + 22}" r="11" class="badge"/>'
        f'<text x="{LX + 22}" y="{y + 26}" class="badge-text">{i + 1}</text>'
    )
    out.append(text(LX + 42, y + 27, f"biome-finder {cmd}", "title"))
    out.append(text(LX + LW - 14, y + 26, status, "status " + status, ' text-anchor="end"'))
    ty = y + 52
    for d in does:
        out.append(text(LX + 22, ty, d, "detail"))
        ty += LINE
    out.append(text(LX + 22, ty + 10, "→ " + output, "output"))
    out.append(f'<rect x="{RX}" y="{y}" width="{RW}" height="{h}" rx="6" class="sources"/>')
    out.append(
        f'<line x1="{RX}" y1="{y + h / 2}" x2="{LX + LW + 4}" y2="{y + h / 2}" class="flow" '
        f'marker-end="url(#arrow)"/>'
    )
    ty = y + 26
    for r in reads:
        out.append(text(RX + 16, ty, "• " + r, "source"))
        ty += LINE
    if skips:
        ty += 6
        for s in skips:
            out.append(text(RX + 16, ty, s, "skip"))
            ty += LINE
    y += h
    cx = LX + LW // 2
    label = HUMAN[i] if i < len(HUMAN) else None
    out.append(f'<line x1="{cx}" y1="{y}" x2="{cx}" y2="{y + GAP - 4}" class="flow" marker-end="url(#arrow)"/>')
    if label:
        out.append(f'<rect x="{cx - 250}" y="{y + 18}" width="500" height="26" rx="13" class="human"/>')
        out.append(text(cx, y + 36, label, "human-text", ' text-anchor="middle"'))
    y += GAP
out.append(f'<rect x="{LX}" y="{y}" width="{RX + RW - LX}" height="44" rx="6" fill="#1f2933"/>')
out.append(text((LX + RX + RW) // 2, y + 27, DOWNSTREAM, "downstream", ' text-anchor="middle"'))
height = y + 44 + 20

svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="{height}" viewBox="0 0 1100 {height}" font-family="Helvetica, Arial, sans-serif">
  <title>biome-finder steps and the data each one reads</title>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#7a8a99"/></marker>
    <style>
      .box {{ fill: #f5f8fb; stroke: #c9d6e3; stroke-width: 1.2; }}
      .sources {{ fill: #fbfbf8; stroke: #d8d8cf; stroke-width: 1.2; }}
      .human {{ fill: #fff7e6; stroke: #E69F00; stroke-dasharray: 4 3; }}
      .badge {{ fill: #0072B2; }}
      .badge-text {{ fill: #ffffff; font-size: 11px; font-weight: bold; text-anchor: middle; }}
      .title {{ font-size: 15px; font-weight: bold; fill: #1f2933; }}
      .heading {{ font-size: 14px; font-weight: bold; fill: #1f2933; }}
      .status {{ font-size: 11.5px; font-weight: bold; }}
      .implemented {{ fill: #009E73; }}
      .planned {{ fill: #8a94a0; }}
      .detail {{ font-size: 12.5px; fill: #44515e; }}
      .output {{ font-family: Menlo, Consolas, "DejaVu Sans Mono", monospace; font-size: 12px; fill: #0b4f8a; }}
      .source {{ font-size: 12px; fill: #1f2933; }}
      .skip {{ font-size: 12px; fill: #D55E00; font-weight: bold; }}
      .human-text {{ font-size: 12px; fill: #7a5200; font-style: italic; }}
      .downstream {{ font-size: 13px; fill: #ffffff; font-weight: bold; }}
      .flow {{ stroke: #7a8a99; stroke-width: 1.5; }}
    </style>
  </defs>
  <rect width="1100" height="{height}" fill="#ffffff"/>
  <text x="{LX}" y="40" class="heading">steps</text>
  <text x="{RX}" y="40" class="heading">data it reads</text>
  {chr(10).join("  " + l for l in out).strip()}
</svg>
"""
open(sys.argv[1], "w").write(svg)

"""Sphinx configuration for biome-finder."""

project = "biome-finder"
author = "biome-finder contributors"
release = "0.1.0"

extensions = ["myst_parser"]
myst_heading_anchors = 4  # links to section headings, e.g. survey.md#worked-example-marine-water-euk
exclude_patterns = ["_build"]
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]  # table cells wrap instead of scrolling sideways

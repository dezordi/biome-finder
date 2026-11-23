"""Clients for the third-party services biome-finder calls, one module per service.

Each client takes a query model (or a plain identifier) from `biome_finder.models` and returns
`Response[T]`. Command modules (`survey`, `search`, …) reach these services only through here."""

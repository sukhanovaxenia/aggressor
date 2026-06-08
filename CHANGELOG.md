# Changelog

All notable changes to AGGRESSOR are documented here.
This project adheres to [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-06-08

### Added
- Installable package with `src/aggressor/` multilevel layout
  (`core`, `rules`, `analysis`, `mutagenesis`, `io`, `cli`).
- `pyproject.toml` (PEP 621), console entry point `aggressor`, `python -m aggressor`.
- `analysis/gatekeepers.py`: automatic, geometry-aware gatekeeper-slot
  detection with optional `--max-gatekeepers-per-apr` budget cap; slots
  ranked by predicted aggregation-propensity reduction; already-gatekeeper
  positions skipped.
- `io/fasta.py`: single source of truth for FASTA I/O and region parsing.
- Test suite (regression + determinism + gatekeeper unit tests), CI, docs.

### Fixed
- **Deterministic output.** Replaced `list(set(...))` mutation ordering and
  the final sort key so identical inputs yield byte-identical output across
  runs (previously varied with the process hash seed).
- **Restored `nearby_hydrophobic_distances`** on `HydrophobicAromaticCluster`
  (dropped during the original module split), with distance + residue
  re-added to its serialization — recovers van der Waals contact-feasibility
  annotation for the 1-pair-plus-hydrophobic condition.
- Removed the duplicated `parse_region` / `_parse_region` (CLI vs engine).

### Changed
- Removed dead `--min-agg-score` argument (parsed but never used).

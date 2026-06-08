# AGGRESSOR

Rule-based *in silico* mutagenesis targeting aggregation-prone regions (APRs).

AGGRESSOR identifies APRs with four physicochemical rules, merges overlapping
clusters into multi-rule hotspots, and proposes substitutions classified by
structural context — with automatic, geometry-aware gatekeeper placement.

- **Deterministic**: identical input yields byte-identical output.
- **No runtime dependencies**: pure standard library, Python >= 3.10.
- **Library + CLI**: `aggressor ...`, `python -m aggressor ...`, or `import aggressor`.

See [Usage](usage.md) to get started and [Biology & methods](biology.md)
for the scientific rationale.

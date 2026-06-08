# AGGRESSOR

**A**ggregation-**G**uided **G**eneration of **RE**gion-specific **S**ubstitution-**OR**iented mutations.

AGGRESSOR performs rule-based *in silico* mutagenesis of protein sequences,
targeting aggregation-prone regions (APRs) and proposing biologically
grounded substitutions — including automatic, geometry-aware **gatekeeper**
placement at APR boundaries.

It implements four physicochemical APR rules (hydrophobic-aliphatic,
aromatic, amide, hydrophobic-aromatic adjacency), merges overlapping
clusters via Union-Find into multi-rule convergence hotspots, classifies
every mutation by structural context (BETA_CORE / GATEKEEPER / BOUNDARY /
FLANKING / DIRECT / INSERTION), and generates single and multi-point
variants.

> References: Rousseau et al., *J Mol Biol* 2006 (gatekeeper hypothesis);
> Beerten et al., *FEBS Lett* 2012 (APR boundary effects);
> Tartaglia et al., *J Mol Biol* 2008 (aggregation propensity scale).

## Installation

```bash
# From source (editable, for development)
pip install -e .[dev]

# From a checkout / wheel
pip install .

# From git
pip install "git+https://github.com/your-org/aggressor.git"

# From bioconda (once published)
conda install -c bioconda -c conda-forge aggressor-mutagenesis
```

Pure standard library — no runtime dependencies. Python >= 3.10.

## Usage

```bash
# Rule-based mutagenesis across a region
aggressor protein.fasta --regions 55:135

# Aggregation analysis only (no mutations)
aggressor protein.fasta --regions all --agg-only

# Restrict to specific rules
aggressor protein.fasta --regions 10:30 --rules hydrophobic_aliphatic aromatic

# Double and triple mutations, parallelised
aggressor protein.fasta --regions 10:30 --multi-mutations 2 3 --threads 4

# Cap gatekeepers per APR for a synthesis budget (otherwise automatic)
aggressor protein.fasta --regions 55:135 --max-gatekeepers-per-apr 1
```

Also available as a module and a library:

```bash
python -m aggressor protein.fasta --regions 55:135
```

```python
from aggressor import mutate_sequence, analyze_region

analysis = analyze_region(seq, start=55, stop=135)
mutations, region_analyses = mutate_sequence(
    seq, positions=[], mutations=["P", "G", "D", "K"], regions=["55:135"],
)
```

## Gatekeeper design philosophy

The number of effective gatekeeper slots is a **structural property of the
APR** (its flanks), not a free integer. AGGRESSOR therefore selects gatekeeper
positions automatically from APR geometry by default. `--max-gatekeepers-per-apr`
is an optional budget cap only: when set, boundary slots are ranked by their
predicted reduction in intrinsic aggregation propensity (Tartaglia scale) and
slots already occupied by a gatekeeper residue are skipped.

## Package layout

```
src/aggressor/
├── core/         # config constants + data models (Cluster, MutationType, ...)
├── rules/        # rule evaluators + registry
├── analysis/     # Union-Find clustering, region analysis, gatekeeper selection
├── mutagenesis/  # single- and multi-point mutation engine
├── io/           # FASTA read/write + region parsing (single source of truth)
└── cli/          # argument parser + application entry point
```

## Development

```bash
pip install -e .[dev]
pytest                 # run the test suite
ruff check src tests   # lint
mkdocs serve           # preview docs at http://127.0.0.1:8000
```

## License

MIT — see `LICENSE`.

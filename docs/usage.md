# Usage

## Command line

```bash
aggressor protein.fasta --regions 55:135
aggressor protein.fasta --regions all --agg-only
aggressor protein.fasta --regions 10:30 --multi-mutations 2 3 --threads 4
aggressor protein.fasta --regions 55:135 --max-gatekeepers-per-apr 1
```

## Library

```python
from aggressor import analyze_region, mutate_sequence

analysis = analyze_region(seq, start=55, stop=135)
mutations, _ = mutate_sequence(
    seq, positions=[], mutations=["P", "G", "D", "K"], regions=["55:135"],
)
```

## Output

Single mutations are written to one FASTA file; multi-point mutations are
organised under `<output>/<level>_mutations/` and split into
`single_region`, `multi_region`, `all_gatekeeper`, `all_core`, and `mixed`.

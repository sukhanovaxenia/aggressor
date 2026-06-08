# Biology & methods

## Aggregation-prone regions

AGGRESSOR detects APRs with four rules, each grounded in the physicochemistry
of cross-beta assembly:

| Rule | Residues | Trigger | Score |
|------|----------|---------|-------|
| hydrophobic_aliphatic | V I L A M | >=3 within 4 positions | 3 |
| aromatic | F Y W | >=2 within 3 positions | 2 |
| amide | Q N | >=2 within 3 positions | 1 |
| hydrophobic_and_aromatic | V I L A M / F Y W | >=2 adjacent pairs, or 1 pair + nearby hydrophobic | 2 |

Overlapping clusters from different rules are merged (Union-Find) into
**multi-rule clusters** whose combined score reflects compounded risk.

## Gatekeepers

Natural proteins suppress aggregation by enriching gatekeeper residues at APR
**flanks**: charged residues (D, E, K, R) impose electrostatic repulsion
between in-register stacking molecules, and proline breaks beta-strand
propagation (Rousseau et al. 2006; Beerten et al. 2012).

A consequence for design: the number of effective gatekeeper slots is a
*structural property* of the APR (its two flanks), not a free integer.
AGGRESSOR therefore chooses gatekeeper positions automatically. The optional
`--max-gatekeepers-per-apr` cap exists only for synthesis-budget reasons;
when set, slots are ranked by predicted reduction in intrinsic aggregation
propensity (Tartaglia et al. 2008) and positions already occupied by a
gatekeeper residue are skipped.

## References

- Rousseau, Schymkowitz & Serrano. *J Mol Biol* 2006.
- Beerten et al. *FEBS Lett* 2012.
- Tartaglia et al. *J Mol Biol* 2008.

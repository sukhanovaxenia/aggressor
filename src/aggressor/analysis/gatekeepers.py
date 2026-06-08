"""
Gatekeeper-slot detection and optimal selection.

Biological rationale
--------------------
Aggregation-prone regions (APRs) form the cross-beta core of amyloid
fibrils. Natural proteins suppress this by enriching *gatekeeper* residues
at the APR flanks (Rousseau et al., J Mol Biol 2006; Beerten et al.,
FEBS Lett 2012; Reumers/Maurer-Stroh and colleagues on sequence-encoded
aggregation safeguards). Two mechanisms dominate:

  * Charged residues (D, E, K, R) impose electrostatic repulsion between
    molecules stacking in-register in the fibril, and raise local
    solubility.
  * Proline is a beta-strand breaker: its backbone cannot donate an amide
    hydrogen bond and its phi is conformationally restricted, so it
    truncates beta-sheet propagation.

A key consequence for *design*: the number of effective gatekeeper slots is
a structural property of the APR — essentially its two boundaries (and the
immediately flanking positions) — not a free integer the user picks. This
module therefore makes the *automatic* choice the primary mode:

  1. Detect whether a gatekeeper substitution at a boundary is even useful
     (skip positions already occupied by a gatekeeper residue — "is it
     possible / needed").
  2. When the caller supplies a budget cap, rank candidate slots by the
     predicted reduction in intrinsic aggregation propensity (Tartaglia et
     al., J Mol Biol 2008 scale) and keep the highest-value ones.

Public API
----------
select_gatekeeper_positions(...)  -> Set[int]   # boundary slots to gatekeep
flanking_gatekeeper_positions(...) -> List[int] # conservative out-of-APR flanks
rank_gatekeeper_slots(...)        -> List[(pos, gain)]
best_gatekeeper_for_position(...) -> Optional[str]
"""
from typing import List, Set, Optional, Tuple, Iterable

from aggressor.core.config import (
    AGGREGATION_PROPENSITY,
    CANONICAL_GATEKEEPER_AAS,
    logger,
)


def _is_gatekeeper_residue(aa: str, gatekeeping_aas: Iterable[str]) -> bool:
    """True if `aa` already acts as a gatekeeper (no benefit to re-introducing)."""
    return aa.upper() in {g.upper() for g in gatekeeping_aas}


def best_gatekeeper_for_position(
        position: int,
        sequence: str,
        gatekeeping_aas: Iterable[str],
) -> Optional[str]:
    """
    Choose the single most propensity-lowering gatekeeper AA for a position.

    Returns the gatekeeper amino acid that yields the lowest intrinsic
    aggregation propensity at this site, or None if the original residue is
    already a gatekeeper (substitution would not help).

    Example:
        >>> best_gatekeeper_for_position(1, "I", {"P", "K", "R", "D", "E"})
        'K'
    """
    original = sequence[position - 1].upper()
    if _is_gatekeeper_residue(original, gatekeeping_aas):
        return None
    candidates = [g.upper() for g in gatekeeping_aas if g.upper() != original]
    if not candidates:
        return None
    # Lowest propensity wins (most aggregation-suppressing)
    return min(candidates, key=lambda g: AGGREGATION_PROPENSITY.get(g, 0.0))


def rank_gatekeeper_slots(
        candidate_positions: List[int],
        sequence: str,
        gatekeeping_aas: Iterable[str],
) -> List[Tuple[int, float]]:
    """
    Rank candidate gatekeeper positions by predicted propensity reduction.

    The "gain" at a position is (propensity of the original residue) minus
    (propensity of the best gatekeeper substitution). High-propensity
    boundary residues (I, V, L, F) yield the largest gains and are therefore
    the most valuable gatekeeper slots.

    Positions whose residue is already a gatekeeper are excluded.

    Returns:
        List of (position, gain) sorted by gain descending, then position
        ascending (deterministic tie-break).
    """
    scored: List[Tuple[int, float]] = []
    for pos in candidate_positions:
        original = sequence[pos - 1].upper()
        gk = best_gatekeeper_for_position(pos, sequence, gatekeeping_aas)
        if gk is None:
            continue
        gain = AGGREGATION_PROPENSITY.get(original, 0.0) - AGGREGATION_PROPENSITY.get(gk, 0.0)
        scored.append((pos, gain))
    scored.sort(key=lambda t: (-t[1], t[0]))
    return scored


def select_gatekeeper_positions(
        cluster_positions: List[int],
        sequence: str,
        gatekeeping_aas: Iterable[str] = CANONICAL_GATEKEEPER_AAS,
        max_gatekeepers: Optional[int] = None,
) -> Set[int]:
    """
    Decide which APR boundary positions should receive gatekeeping AAs.

    Default (``max_gatekeepers is None``): both cluster boundaries — i.e.
    {min, max} — exactly the automatic behaviour. This is intentional: the
    two flanks are the biologically meaningful gatekeeper slots, so the
    count is detected from APR geometry, not requested by the user.

    With a cap: the boundary slots are ranked by predicted propensity
    reduction and only the top ``max_gatekeepers`` are returned. Slots whose
    residue is already a gatekeeper are skipped, so the cap is applied to
    *useful* slots only.

    Args:
        cluster_positions: All positions in the merged APR cluster.
        sequence: Full protein sequence (for residue lookup).
        gatekeeping_aas: Allowed gatekeeper amino acids.
        max_gatekeepers: Optional budget cap on gatekept slots per APR.

    Returns:
        Set of 1-indexed positions that should receive gatekeeping AAs.

    Example:
        >>> sorted(select_gatekeeper_positions([5, 6, 7], "AAAAIFLAA"))
        [5, 7]
    """
    if not cluster_positions:
        return set()

    cluster_min = min(cluster_positions)
    cluster_max = max(cluster_positions)
    boundary_slots = {cluster_min, cluster_max}

    if max_gatekeepers is None:
        # Automatic mode: geometry decides the count.
        return boundary_slots

    if max_gatekeepers <= 0:
        return set()

    ranked = rank_gatekeeper_slots(sorted(boundary_slots), sequence, gatekeeping_aas)
    chosen = {pos for pos, _ in ranked[:max_gatekeepers]}
    if not chosen:
        logger.debug(
            "No useful gatekeeper slot at APR %s:%s "
            "(boundaries already gatekeeper residues)",
            cluster_min, cluster_max,
        )
    return chosen


def flanking_gatekeeper_positions(
        cluster_positions: List[int],
        region_start: int,
        region_end: int,
) -> List[int]:
    """
    Return the positions immediately flanking an APR (one N-terminal, one
    C-terminal), clamped to the analysis region.

    This supports the conservative gatekeeper-design strategy in which a
    charge is introduced *just outside* the APR core, leaving the cross-beta
    register intact while walling it off electrostatically — distinct from
    substituting a terminal core residue. Exposed for downstream/expert use;
    not part of the default mutation path so that standard output stays
    backward-compatible.
    """
    if not cluster_positions:
        return []
    flanks = []
    left = min(cluster_positions) - 1
    right = max(cluster_positions) + 1
    if left >= region_start:
        flanks.append(left)
    if right <= region_end:
        flanks.append(right)
    return flanks

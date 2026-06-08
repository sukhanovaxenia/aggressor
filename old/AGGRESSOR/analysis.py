"""
Region analysis for aggregation-prone cluster detection.

Provides the main analysis pipeline that:
1. Applies multiple aggregation rules to a sequence region
2. Merges overlapping clusters using Union-Find
3. Identifies multi-rule convergence hotspots
4. Generates structured analysis reports

This module orchestrates the clustering rules and merging algorithms
to produce a comprehensive aggregation propensity profile for
any protein region.

Architecture:
    analyze_region() — main entry point
    classify_position_context() — position classification relative to clusters
    classify_mutation_type() — mutation type determination

Dependencies:
    - rules.py: Rule evaluators and registry
    - clustering.py: Cluster merging algorithms
    - models.py: Data structures
    - config.py: Constants and thresholds
"""
from typing import Dict, Any, Optional, List, Set

from config import (
    MAX_GAP_FOR_MERGING,
    GATEKEEPER_DISTANCE,
    FLANKING_DISTANCE,
    CANONICAL_GATEKEEPER_AAS,
    logger
)
from models import (
    Cluster,
    MultiRuleCluster,
    MutationType,
    PositionContext,
    ClusterType
)
from rules import RuleRegistry, DEFAULT_REGISTRY
from clustering import merge_overlapping_clusters_unionfind


# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

def analyze_region(
        sequence: str,
        start: int,
        stop: int,
        registry: RuleRegistry = None,
        selected_rules: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Analyze a protein region for aggregation-prone clusters.

    This is the main analysis function that orchestrates the entire
    aggregation detection pipeline:
    1. Apply all (or selected) rules to the region
    2. Collect identified clusters
    3. Merge overlapping clusters using Union-Find
    4. Identify multi-rule convergence hotspots
    5. Build comprehensive results structure

    Args:
        sequence: Full protein sequence (1-letter codes)
        start: Region start position (1-indexed, inclusive)
        stop: Region end position (1-indexed, inclusive)
        registry: RuleRegistry with evaluators to use.
                 If None, uses DEFAULT_REGISTRY.
        selected_rules: Optional list of rule names to apply.
                       If None, all registered rules are used.

    Returns:
        Dictionary with comprehensive analysis results:
        {
            'region': (start, stop),
            'sequence': region_sequence,
            'length': region_length,
            'rules': {rule_name: rule_results, ...},
            'merged_clusters': [cluster_dicts, ...],
            'multi_rule_clusters': [cluster_dicts, ...],
            'aggregation_hotspots': [positions, ...],
            'all_cluster_positions': [positions, ...]
        }

    Example:
        >>> seq = "ALVIFQNWYG"
        >>> results = analyze_region(seq, 1, 10)
        >>> results['region']
        (1, 10)
        >>> len(results['merged_clusters'])
        2
    """
    if registry is None:
        registry = DEFAULT_REGISTRY

    region_seq = sequence[start - 1:stop]

    logger.debug(f"Analyzing region {start}:{stop} ({len(region_seq)} residues)")
    logger.debug(f"Region sequence: {region_seq}")

    # Step 1: Evaluate all applicable rules
    rule_results = registry.evaluate_region(
        sequence, start, stop, selected_rules
    )

    # Step 2: Collect all clusters across rules
    all_clusters: List[Cluster] = []
    for rule_name, clusters in rule_results.items():
        all_clusters.extend(clusters)

    logger.debug(f"Total raw clusters found: {len(all_clusters)}")

    # Step 3: Merge overlapping clusters using Union-Find
    merged_clusters = merge_overlapping_clusters_unionfind(
        all_clusters, sequence, gap_tolerance=MAX_GAP_FOR_MERGING
    )

    logger.debug(f"Clusters after merging: {len(merged_clusters)}")

    # Step 4: Identify multi-rule convergence clusters
    multi_rule_clusters = [
        c for c in merged_clusters
        if isinstance(c, MultiRuleCluster)
    ]

    if multi_rule_clusters:
        logger.debug(
            f"Multi-rule convergence clusters: {len(multi_rule_clusters)}"
        )
        for mrc in multi_rule_clusters:
            logger.debug(
                f"  Rules: {mrc.rules}, "
                f"Score: {mrc.combined_aggregation_score}"
            )

    # Step 5: Collect all hotspot positions
    hotspot_positions: Set[int] = set()
    for cluster in merged_clusters:
        hotspot_positions.update(cluster.positions)

    # Step 6: Build comprehensive results structure
    results = {
        'region': (start, stop),
        'sequence': region_seq,
        'length': len(region_seq),
        'rules': {},
        'merged_clusters': [c.to_dict() for c in merged_clusters],
        'multi_rule_clusters': [c.to_dict() for c in multi_rule_clusters],
        'aggregation_hotspots': sorted(hotspot_positions),
        'all_cluster_positions': sorted(hotspot_positions),
    }

    # Step 7: Add per-rule details with backward-compatible fields
    for rule_name, clusters in rule_results.items():
        evaluator = registry.get(rule_name)

        rule_entry = {
            'description': evaluator.description,
            'qualifying_clusters': [c.to_dict() for c in clusters],
            'condition_met': len(clusters) > 0,
            'aggregation_score': evaluator.aggregation_score,
            # Backward compatibility fields
            'matching_positions': [],
            'matching_residues': [],
            'clusters': [],
        }

        # Populate backward-compatible fields
        for cluster in clusters:
            rule_entry['matching_positions'].extend(cluster.positions)
            rule_entry['matching_residues'].extend(cluster.residues)
            rule_entry['clusters'].append(list(cluster.positions))

        # Remove duplicates while preserving order
        rule_entry['matching_positions'] = list(
            dict.fromkeys(rule_entry['matching_positions'])
        )
        rule_entry['matching_residues'] = list(
            dict.fromkeys(rule_entry['matching_residues'])
        )

        # Special handling for hydrophobic_and_aromatic rule
        if rule_name == 'hydrophobic_and_aromatic' and clusters:
            rule_entry['special_clusters'] = [
                c.to_dict() for c in clusters
            ]

        results['rules'][rule_name] = rule_entry

    return results


# =============================================================================
# POSITION CLASSIFICATION
# =============================================================================

def classify_position_context(
        position: int,
        clusters: List[Dict[str, Any]],
        gatekeeper_distance: int = GATEKEEPER_DISTANCE,
        flanking_distance: int = FLANKING_DISTANCE
) -> PositionContext:
    """
    Classify a position's structural context relative to aggregation clusters.

    Determines whether a given position falls within, at the boundary of,
    or near an aggregation-prone region (APR). This classification guides
    mutation strategy selection.

    Classification hierarchy:
    - BETA_CORE: Position is within an identified APR cluster
    - GATEKEEPER: Within gatekeeper_distance of an APR boundary
    - FLANKING: Within flanking_distance but beyond gatekeeper zone
    - UNKNOWN: Beyond all thresholds or no clusters present

    Args:
        position: 1-indexed sequence position to classify
        clusters: List of cluster dictionaries from analyze_region()
        gatekeeper_distance: Max distance from APR boundary for gatekeeper
        flanking_distance: Max distance from APR for flanking region

    Returns:
        PositionContext with classification and distance metrics

    Example:
        >>> clusters = [{'positions': [5, 6, 7]}]
        >>> ctx = classify_position_context(4, clusters)
        >>> ctx.context_type == MutationType.GATEKEEPER
        True
        >>> ctx.distance_to_nearest_apr
        1
    """
    if not clusters:
        return PositionContext(
            position=position,
            context_type=MutationType.UNKNOWN,
            distance_to_nearest_apr=999999,
        )

    min_distance = 999999
    nearest_boundary = None
    containing_cluster_id = None

    for i, cluster in enumerate(clusters):
        positions = cluster.get('positions', [])
        if not positions:
            continue

        cluster_min = min(positions)
        cluster_max = max(positions)

        # Check if position is within the cluster
        if position in positions:
            return PositionContext(
                position=position,
                context_type=MutationType.BETA_CORE,
                distance_to_nearest_apr=0,
                nearest_apr_boundary=position,
                apr_id=i,
            )

        # Calculate distance to nearest cluster boundary
        if position < cluster_min:
            distance = cluster_min - position
            boundary = cluster_min
        elif position > cluster_max:
            distance = position - cluster_max
            boundary = cluster_max
        else:
            # Position is between cluster positions but not IN cluster
            # Find the closest position
            distance = min(abs(position - p) for p in positions)
            boundary = min(positions, key=lambda p: abs(position - p))

        if distance < min_distance:
            min_distance = distance
            nearest_boundary = boundary
            containing_cluster_id = i

    # Classify based on distance thresholds
    if min_distance <= gatekeeper_distance:
        context_type = MutationType.GATEKEEPER
    elif min_distance <= flanking_distance:
        context_type = MutationType.FLANKING
    else:
        context_type = MutationType.UNKNOWN

    return PositionContext(
        position=position,
        context_type=context_type,
        distance_to_nearest_apr=min_distance,
        nearest_apr_boundary=nearest_boundary,
        apr_id=(
            containing_cluster_id
            if min_distance <= flanking_distance
            else None
        ),
    )


def classify_mutation_type(
        position: int,
        new_aa: str,
        cluster_positions: List[int],
        cluster_min: int,
        cluster_max: int,
        gatekeeper_aas: Optional[Set[str]] = None
) -> MutationType:
    """
    Determine mutation type based on position context AND amino acid properties.

    A true gatekeeper mutation requires BOTH:
    1. Position at or near APR boundary (within GATEKEEPER_DISTANCE)
    2. Mutation TO a canonical gatekeeper amino acid (P, K, R, D, E)

    This dual requirement ensures biological accuracy: only mutations
    that can actually disrupt β-sheet extension at cluster boundaries
    are classified as gatekeepers.

    Args:
        position: 1-indexed position of the mutation
        new_aa: Single-letter code of the new amino acid
        cluster_positions: All positions in the cluster
        cluster_min: N-terminal boundary of the cluster
        cluster_max: C-terminal boundary of the cluster
        gatekeeper_aas: Set of amino acids classified as gatekeepers.
                       If None, uses CANONICAL_GATEKEEPER_AAS.

    Returns:
        Appropriate MutationType classification

    Example:
        >>> classify_mutation_type(5, 'P', [5,6,7], 5, 7)
        <MutationType.GATEKEEPER: 2>
        >>> classify_mutation_type(6, 'G', [5,6,7], 5, 7)
        <MutationType.BETA_CORE: 1>
    """
    if gatekeeper_aas is None:
        gatekeeper_aas = CANONICAL_GATEKEEPER_AAS

    is_boundary = (position == cluster_min or position == cluster_max)
    is_adjacent_to_boundary = (
            abs(position - cluster_min) <= GATEKEEPER_DISTANCE or
            abs(position - cluster_max) <= GATEKEEPER_DISTANCE
    )
    is_gatekeeper_aa = new_aa.upper() in gatekeeper_aas

    if position in cluster_positions:
        # Position is within the cluster
        if is_boundary and is_gatekeeper_aa:
            return MutationType.GATEKEEPER
        elif is_boundary:
            return MutationType.BOUNDARY
        else:
            return MutationType.BETA_CORE
    elif is_adjacent_to_boundary:
        if is_gatekeeper_aa:
            return MutationType.GATEKEEPER
        else:
            return MutationType.FLANKING
    else:
        return MutationType.FLANKING
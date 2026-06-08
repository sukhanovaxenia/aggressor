"""
Rule definitions and registry for aggregation-prone region detection.

Implements the Strategy pattern for rule-based identification of
aggregation-prone regions (APRs) in protein sequences. Each rule
detects specific physicochemical patterns associated with amyloid
formation.

Rules detect clusters based on:
- Hydrophobic aliphatic residues (V, I, L, A, M)
- Aromatic residues (F, Y, W) 
- Amide-containing residues (Q, N)
- Hydrophobic-aromatic adjacency patterns

Architecture:
    BaseClusterEvaluator (ABC-like with common clustering logic)
    ├── HydrophobicAliphaticEvaluator
    ├── AromaticEvaluator
    ├── AmideEvaluator
    └── HydrophobicAromaticEvaluator (special pair-detection logic)

    RuleRegistry: Container managing evaluator registration and lookup

References:
- Rousseau et al., J Mol Biol 2006 (gatekeeper hypothesis)
- Alberti et al., Cell 2009 (Q/N-rich prion domains)
- King et al., Science 2012 (amide zipper formation)
- Tartaglia et al., J Mol Biol 2008 (aggregation propensity)
"""
from typing import List, Dict, Optional

from config import logger, DEFAULT_MUTATIONS
from models import Cluster, HydrophobicAromaticCluster, ClusterEvaluator


# =============================================================================
# BASE CLUSTER EVALUATOR
# =============================================================================

class BaseClusterEvaluator:
    """
    Base implementation providing common clustering logic.

    Uses the Template Method pattern: subclasses override class-level
    attributes to define rule-specific behavior while inheriting the
    core clustering algorithm.

    Subclassing Guide:
        1. Set NAME, DESCRIPTION, RESIDUES (class constants)
        2. Override MIN_CLUSTER_SIZE if different from default (2)
        3. Override MAX_GAP if different from default (3)
        4. Override AGGREGATION_SCORE for rule weighting
        5. Override find_clusters() only for special detection logic

    Attributes:
        NAME: Unique rule identifier string
        DESCRIPTION: Human-readable rule description
        RESIDUES: FrozenSet of target amino acid codes
        MIN_CLUSTER_SIZE: Minimum residues required for a cluster
        MAX_GAP: Maximum positions between cluster members
        AGGREGATION_SCORE: Integer weight for this rule
    """
    # Override these in subclasses
    NAME: str = "base"
    DESCRIPTION: str = "Base evaluator"
    RESIDUES: frozenset = frozenset()
    MIN_CLUSTER_SIZE: int = 2
    MAX_GAP: int = 3
    AGGREGATION_SCORE: int = 1

    @property
    def name(self) -> str:
        """Unique identifier for this rule evaluator."""
        return self.NAME

    @property
    def description(self) -> str:
        """Human-readable description of what this rule detects."""
        return self.DESCRIPTION

    @property
    def aggregation_score(self) -> int:
        """Integer weight indicating aggregation contribution."""
        return self.AGGREGATION_SCORE

    @property
    def residues(self) -> frozenset:
        """Set of amino acid single-letter codes this rule operates on."""
        return self.RESIDUES

    def find_clusters(
            self,
            sequence: str,
            start: int,
            stop: int
    ) -> List[Cluster]:
        """
        Standard clustering implementation using sliding window.

        Algorithm:
        1. Extract region from full sequence
        2. Find all positions matching target residues
        3. Group positions within MAX_GAP into clusters
        4. Filter clusters by MIN_CLUSTER_SIZE
        5. Return Cluster objects with absolute positions

        Args:
            sequence: Full protein sequence
            start: 1-indexed start of region to analyze
            stop: 1-indexed end of region to analyze

        Returns:
            List of Cluster objects meeting criteria
        """
        region_seq = sequence[start - 1:stop]

        # Find positions of all matching residues in the region
        matching_positions = [
            start + i
            for i, aa in enumerate(region_seq)
            if aa in self.RESIDUES
        ]

        # Group into clusters based on proximity
        raw_clusters = self._cluster_positions(matching_positions)

        # Create Cluster objects for qualifying clusters
        return [
            Cluster(
                positions=tuple(cluster),
                residues=tuple(sequence[p - 1] for p in cluster),
                rule_name=self.NAME,
                aggregation_score=self.AGGREGATION_SCORE
            )
            for cluster in raw_clusters
            if len(cluster) >= self.MIN_CLUSTER_SIZE
        ]

    def _cluster_positions(self, positions: List[int]) -> List[List[int]]:
        """
        Group positions into clusters based on MAX_GAP threshold.

        Consecutive positions within MAX_GAP are grouped together.
        This is a single-pass greedy clustering algorithm.

        Args:
            positions: Sorted list of 1-indexed positions

        Returns:
            List of clusters, each being a list of positions
        """
        if not positions:
            return []

        clusters = []
        current = [positions[0]]

        for pos in positions[1:]:
            if pos - current[-1] <= self.MAX_GAP:
                # Extend current cluster
                current.append(pos)
            else:
                # Save current if it meets minimum size
                if len(current) >= 2:
                    clusters.append(current)
                # Start new cluster
                current = [pos]

        # Don't forget the last cluster
        if len(current) >= 2:
            clusters.append(current)

        return clusters


# =============================================================================
# CONCRETE RULE EVALUATORS
# =============================================================================

class HydrophobicAliphaticEvaluator(BaseClusterEvaluator):
    """
    Detects clusters of hydrophobic aliphatic residues (V, I, L, A, M).

    Biological basis:
        Aliphatic side chains drive aggregation through hydrophobic
        collapse and subsequent β-sheet formation. The β-branched
        residues (V, I) have particularly high intrinsic aggregation
        propensity due to their conformational preferences favoring
        extended (β-strand) structures.

    Triggering criteria:
        ≥3 hydrophobic aliphatic residues within a span of 4 positions.
        This strict threshold reduces false positives while capturing
        genuine aggregation-prone segments.

    Aggregation score: 3 (highest single-rule score)
    """
    NAME = "hydrophobic_aliphatic"
    DESCRIPTION = "V, I, L, A, M clustered (≥3 residues within 4 positions)"
    RESIDUES = frozenset('VILAM')
    MIN_CLUSTER_SIZE = 3
    MAX_GAP = 4
    AGGREGATION_SCORE = 3


class AromaticEvaluator(BaseClusterEvaluator):
    """
    Detects clusters of aromatic residues (F, Y, W).

    Biological basis:
        Aromatic residues contribute to aggregation through π-π stacking
        interactions and hydrophobic contacts. Phenylalanine is particularly
        aggregation-prone due to its high hydrophobicity. Tyrosine can
        modulate aggregation through its hydroxyl group's hydrogen bonding
        capability. Tryptophan's large indole ring creates strong
        hydrophobic and π-interactions.

    Triggering criteria:
        ≥2 aromatic residues within a span of 3 positions.
        Even two adjacent aromatics can nucleate aggregation through
        π-stacking.

    Aggregation score: 2
    """
    NAME = "aromatic"
    DESCRIPTION = "F, Y, W clustered (≥2 residues within 3 positions)"
    RESIDUES = frozenset('FYW')
    MIN_CLUSTER_SIZE = 2
    MAX_GAP = 3
    AGGREGATION_SCORE = 2


class AmideEvaluator(BaseClusterEvaluator):
    """
    Detects clusters of amide-containing residues (Q, N).

    Biological basis:
        Q/N-rich regions are hallmarks of prion-like domains and
        intrinsically disordered regions prone to liquid-liquid phase
        separation (LLPS). The amide side chains form "polar zipper"
        hydrogen bond networks that can nucleate and stabilize cross-β
        structure. These interactions are particularly important in
        yeast prions (e.g., Sup35) and human proteins associated with
        neurodegenerative diseases (e.g., FUS, hnRNPA1).

    Triggering criteria:
        ≥2 Q or N residues within a span of 3 positions.

    References:
        - Alberti et al., Cell 2009 (Q/N-rich prion domains)
        - King et al., Science 2012 (amide zipper formation)

    Aggregation score: 1
    """
    NAME = "amide"
    DESCRIPTION = "Q, N clustered (≥2 residues within 3 positions)"
    RESIDUES = frozenset('QN')
    MIN_CLUSTER_SIZE = 2
    MAX_GAP = 3
    AGGREGATION_SCORE = 1


class HydrophobicAromaticEvaluator(BaseClusterEvaluator):
    """
    Detects hydrophobic-aromatic adjacency patterns using special logic.

    Unlike other evaluators that use simple positional clustering, this
    evaluator implements two distinct detection mechanisms:

    Condition A: ≥2 hydrophobic-aromatic adjacent pairs
        Multiple adjacent H-A or A-H pairs indicate a pattern of
        alternating hydrophobicity that strongly favors β-sheet
        formation.

    Condition B: 1 pair + additional hydrophobic within 3 positions
        A single pair stabilized by a nearby hydrophobic residue
        can also nucleate aggregation.

    Biological basis:
        Adjacent hydrophobic-aromatic residues create particularly
        stable β-strand segments in amyloid fibrils. The aromatic ring
        intercalates between aliphatic chains, optimizing van der Waals
        packing in the fibril core. This motif is common in many
        amyloidogenic proteins including Aβ, IAPP, and α-synuclein.

    Residue classification:
        Hydrophobic: V, I, L, A, M (aliphatic)
        Aromatic: F, Y, W (ring-containing)

    Aggregation score: 2
    """
    NAME = "hydrophobic_and_aromatic"
    DESCRIPTION = "Hydrophobic (V,I,L,A,M) adjacent to aromatic (F,Y,W)"
    RESIDUES = frozenset('VILAMFYW')
    MIN_CLUSTER_SIZE = 2
    MAX_GAP = 1
    AGGREGATION_SCORE = 2

    # Residue classifications for pair detection
    HYDROPHOBIC = frozenset('VILAM')
    AROMATIC = frozenset('FYW')

    def find_clusters(
            self,
            sequence: str,
            start: int,
            stop: int
    ) -> List[Cluster]:
        """
        Specialized detection using pair-based logic.

        Overrides the base implementation because simple positional
        clustering cannot capture the adjacency pattern requirement.

        Algorithm:
        1. Scan region for all adjacent H-A or A-H pairs
        2. Condition A: ≥2 pairs → cluster all pair positions
        3. Condition B: 1 pair + nearby hydrophobic → expand cluster
        4. Deduplicate clusters by position sets

        Args:
            sequence: Full protein sequence
            start: 1-indexed start of region
            stop: 1-indexed end of region

        Returns:
            List of HydrophobicAromaticCluster objects
        """
        region_seq = sequence[start - 1:stop]

        # Step 1: Find all adjacent hydrophobic-aromatic pairs
        pairs = []
        for i in range(len(region_seq) - 1):
            aa1, aa2 = region_seq[i], region_seq[i + 1]
            pos1, pos2 = start + i, start + i + 1

            is_pair = (
                    (aa1 in self.HYDROPHOBIC and aa2 in self.AROMATIC) or
                    (aa1 in self.AROMATIC and aa2 in self.HYDROPHOBIC)
            )
            if is_pair:
                pairs.append({
                    'positions': (pos1, pos2),
                    'residues': (aa1, aa2),
                    'interaction': f"{aa1}{pos1}-{aa2}{pos2}"
                })

        if not pairs:
            return []

        # Collect all hydrophobic positions in the region (for Condition B)
        hydrophobic_positions = {
            start + i
            for i, aa in enumerate(region_seq)
            if aa in self.HYDROPHOBIC
        }

        clusters = []

        # Condition A: ≥2 pairs → merge all pair positions into one cluster
        if len(pairs) >= 2:
            pair_positions = set()
            pair_interactions = []
            for pair in pairs:
                pair_positions.update(pair['positions'])
                pair_interactions.append(pair['interaction'])

            sorted_positions = tuple(sorted(pair_positions))
            cluster_residues = tuple(sequence[p - 1] for p in sorted_positions)

            hydrophobic_count = sum(
                1 for r in cluster_residues if r in self.HYDROPHOBIC
            )
            aromatic_count = sum(
                1 for r in cluster_residues if r in self.AROMATIC
            )

            clusters.append(HydrophobicAromaticCluster(
                positions=sorted_positions,
                residues=cluster_residues,
                rule_name=self.NAME,
                aggregation_score=self.AGGREGATION_SCORE,
                pair_count=len(pairs),
                condition='at_least_2_pairs',
                hydrophobic_count=hydrophobic_count,
                aromatic_count=aromatic_count,
                pair_interactions=tuple(pair_interactions),
                nearby_hydrophobic_positions=()
            ))

        # Condition B: 1 pair + nearby hydrophobic (within 3 positions)
        for pair in pairs:
            pair_position_set = set(pair['positions'])

            # Find hydrophobic residues within 3 positions of the pair
            nearby = []
            for hp in hydrophobic_positions:
                if hp in pair_position_set:
                    continue
                min_distance = min(abs(hp - p) for p in pair['positions'])
                if min_distance <= 3:
                    nearby.append(hp)

            if nearby:
                all_positions = tuple(
                    sorted(pair_position_set | set(nearby))
                )
                cluster_residues = tuple(
                    sequence[p - 1] for p in all_positions
                )

                hydrophobic_count = sum(
                    1 for r in cluster_residues if r in self.HYDROPHOBIC
                )
                aromatic_count = sum(
                    1 for r in cluster_residues if r in self.AROMATIC
                )

                clusters.append(HydrophobicAromaticCluster(
                    positions=all_positions,
                    residues=cluster_residues,
                    rule_name=self.NAME,
                    aggregation_score=self.AGGREGATION_SCORE,
                    pair_count=1,
                    condition='1_pair_plus_hydrophobic',
                    hydrophobic_count=hydrophobic_count,
                    aromatic_count=aromatic_count,
                    pair_interactions=(pair['interaction'],),
                    nearby_hydrophobic_positions=tuple(nearby)
                ))

        # Deduplicate clusters by position sets
        seen = set()
        unique = []
        for cluster in clusters:
            if cluster.positions not in seen:
                seen.add(cluster.positions)
                unique.append(cluster)

        return unique


# =============================================================================
# RULE REGISTRY
# =============================================================================

class RuleRegistry:
    """
    Centralized registry for managing cluster evaluation rules.

    Implements a simple service locator pattern for evaluator instances.
    Supports dynamic registration of new rules without modifying
    existing code (Open/Closed Principle).

    Usage:
        >>> registry = RuleRegistry()
        >>> registry.register(HydrophobicAliphaticEvaluator())
        >>> evaluator = registry.get("hydrophobic_aliphatic")
        >>> clusters = evaluator.find_clusters(seq, 1, 100)

    Features:
    - Named registration with collision detection
    - Rule listing for UI/dynamic selection
    - Batch evaluation of multiple rules on a region
    """

    def __init__(self):
        """Initialize empty registry."""
        self._evaluators: Dict[str, ClusterEvaluator] = {}

    def register(self, evaluator: ClusterEvaluator) -> None:
        """
        Register an evaluator instance under its name.

        If an evaluator with the same name already exists,
        a warning is logged and the old evaluator is replaced.

        Args:
            evaluator: Instance implementing ClusterEvaluator protocol

        Example:
            >>> registry.register(HydrophobicAliphaticEvaluator())
        """
        if evaluator.name in self._evaluators:
            logger.warning(
                f"Overwriting existing evaluator: {evaluator.name}"
            )
        self._evaluators[evaluator.name] = evaluator
        logger.debug(f"Registered evaluator: {evaluator.name}")

    def get(self, name: str) -> ClusterEvaluator:
        """
        Retrieve an evaluator by its registered name.

        Args:
            name: Rule name (case-sensitive, e.g., 'hydrophobic_aliphatic')

        Returns:
            The registered ClusterEvaluator instance

        Raises:
            KeyError: If the rule name is not registered

        Example:
            >>> evaluator = registry.get("hydrophobic_aliphatic")
        """
        if name not in self._evaluators:
            available = list(self._evaluators.keys())
            raise KeyError(
                f"Unknown rule: '{name}'. "
                f"Available: {available}"
            )
        return self._evaluators[name]

    def list_rules(self) -> List[str]:
        """
        Get list of all registered rule names.

        Returns:
            Sorted list of rule name strings

        Example:
            >>> registry.list_rules()
            ['amide', 'aromatic', 'hydrophobic_aliphatic', 'hydrophobic_and_aromatic']
        """
        return list(self._evaluators.keys())

    def evaluate_region(
            self,
            sequence: str,
            start: int,
            stop: int,
            rules: Optional[List[str]] = None
    ) -> Dict[str, List[Cluster]]:
        """
        Evaluate all (or selected) rules on a sequence region.

        This is the main entry point for region analysis. It applies
        each specified rule to the given region and collects results.

        Args:
            sequence: Full protein sequence
            start: 1-indexed start position of region
            stop: 1-indexed end position of region
            rules: Optional list of rule names to apply.
                   If None, all registered rules are applied.

        Returns:
            Dictionary mapping rule_name → list of Cluster objects

        Example:
            >>> results = registry.evaluate_region(seq, 10, 50)
            >>> for rule_name, clusters in results.items():
            ...     print(f"{rule_name}: {len(clusters)} clusters")
        """
        target_rules = rules if rules is not None else self.list_rules()

        results = {}
        for rule_name in target_rules:
            if rule_name not in self._evaluators:
                logger.warning(f"Skipping unknown rule: {rule_name}")
                continue

            evaluator = self._evaluators[rule_name]
            clusters = evaluator.find_clusters(sequence, start, stop)
            results[rule_name] = clusters

            logger.debug(
                f"Rule '{rule_name}': found {len(clusters)} cluster(s) "
                f"in region {start}:{stop}"
            )

        return results


# =============================================================================
# DEFAULT REGISTRY AND COMPATIBILITY LAYER
# =============================================================================

def create_default_registry() -> RuleRegistry:
    """
    Create a RuleRegistry pre-loaded with all standard aggregation rules.

    Returns:
        RuleRegistry with hydrophobic_aliphatic, aromatic, amide,
        and hydrophobic_and_aromatic evaluators registered.

    Example:
        >>> registry = create_default_registry()
        >>> "hydrophobic_aliphatic" in registry.list_rules()
        True
    """
    registry = RuleRegistry()
    registry.register(HydrophobicAliphaticEvaluator())
    registry.register(AromaticEvaluator())
    registry.register(AmideEvaluator())
    registry.register(HydrophobicAromaticEvaluator())
    return registry


# Module-level default registry instance
DEFAULT_REGISTRY = create_default_registry()


def _generate_rules_dict(registry: RuleRegistry) -> Dict[str, Dict]:
    """
    Generate RULES dictionary from registry for backward compatibility.

    Some parts of the codebase may rely on the old dictionary format
    for rule configuration. This function bridges the gap between
    the new object-oriented registry and the legacy dict format.

    Args:
        registry: A RuleRegistry instance with evaluators

    Returns:
        Dictionary with rule configuration in legacy format
    """
    rules = {}
    for rule_name in registry.list_rules():
        evaluator = registry.get(rule_name)
        rules[rule_name] = {
            'description': evaluator.description,
            'residues': set(evaluator.residues),
            'min_cluster_size': getattr(evaluator, 'MIN_CLUSTER_SIZE', 2),
            'max_gap': getattr(evaluator, 'MAX_GAP', 3),
            'mutations': DEFAULT_MUTATIONS,
            'aggregation_score': evaluator.aggregation_score,
        }
    return rules


# Backward-compatible RULES dictionary
RULES = _generate_rules_dict(DEFAULT_REGISTRY)
#!/usr/bin/env python3
"""
AGGRESSOR: Aggregation-Guided Generation of REgion-Specific Substitution ORiented mutations

In Silico Mutagenesis Script with Rule-Based Mutations and Multi-Point Mutation Support.
Performs point mutations and insertions at specified positions in protein sequences
with rule-based mutagenesis in specified regions.

Features:
- Strategy Pattern for rule evaluation
- Union-Find for efficient cluster merging
- Biologically accurate gatekeeper classification
- Multi-point mutation generation with parallelization
- Structured mutation type classification

Usage: python AGGRESSOR.py <input_file> [options]

References:
- Rousseau et al., J Mol Biol 2006 (gatekeeper hypothesis)
- Beerten et al., FEBS Lett 2012 (APR boundary effects)
- Tartaglia et al., J Mol Biol 2008 (aggregation propensity scale)
"""
import argparse
import re
import sys
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum, auto
from multiprocessing import cpu_count
from typing import (
    List, Tuple, Dict, Set, Any, Union, Optional,
    FrozenSet, Protocol, runtime_checkable
)
from pathlib import Path
from itertools import combinations, product

# =============================================================================
# CONSTANTS
# =============================================================================

VALID_AAS = set('ACDEFGHIKLMNPQRSTVWY')
FASTA_LINE_LENGTH = 60
MAX_GAP_FOR_MERGING = 2
MAX_MUTATION_LEVEL = 10  # Hard upper limit for multi-mutations

# Default mutation list - customizable here
DEFAULT_MUTATIONS = ['P', 'G', 'D', 'K']

# Gatekeeping amino acids - charged residues and proline that disrupt β-extension
# These are only applied to positions at cluster boundaries
GATEKEEPING_AAS = ['Y']

# Canonical gatekeeper amino acids (introduce charge or conformational constraint)
# Based on Rousseau et al., J Mol Biol 2006
CANONICAL_GATEKEEPER_AAS = frozenset({'P', 'K', 'R', 'D', 'E'})

# Distance thresholds for position classification
GATEKEEPER_DISTANCE = 3  # Max distance from APR boundary for gatekeeper
FLANKING_DISTANCE = 6    # Max distance for flanking region

# Empirically-derived aggregation propensity values
# Source: Tartaglia et al., J Mol Biol 2008
AGGREGATION_PROPENSITY: Dict[str, float] = {
    'I': 1.822,  'V': 1.594,  'L': 1.380,  'F': 1.376,
    'Y': 0.888,  'W': 0.893,  'M': 0.739,  'A': 0.411,
    'C': 0.382,  'T': 0.039,  'S': -0.228, 'G': -0.535,
    'N': -0.547, 'Q': -0.691, 'H': -0.731, 'P': -1.402,
    'D': -1.836, 'E': -1.892, 'K': -2.030, 'R': -1.814,
}

# Module-level logger
logger = logging.getLogger(__name__)


# =============================================================================
# ENUMS AND TYPE DEFINITIONS
# =============================================================================

class MutationType(Enum):
    """
    Classification of mutations by structural/functional context.
    
    Based on the gatekeeper hypothesis (Rousseau et al., 2006) and
    aggregation-prone region (APR) architecture.
    
    Categories:
        BETA_CORE: Within identified aggregation cluster (highest aggregation risk)
        GATEKEEPER: At cluster boundary with gatekeeper AA (disrupts β-extension)
        FLANKING: Adjacent to APR but not gatekeeper position/AA
        BOUNDARY: At APR boundary but mutation is not to gatekeeper AA
        DIRECT: User-specified position (no rule context)
        INSERTION: Amino acid insertion
        IDR_PROXIMAL: Near intrinsically disordered region (future use)
        UNKNOWN: Cannot classify
    """
    BETA_CORE = auto()
    GATEKEEPER = auto()
    FLANKING = auto()
    BOUNDARY = auto()
    DIRECT = auto()
    INSERTION = auto()
    IDR_PROXIMAL = auto()
    UNKNOWN = auto()

    @classmethod
    def from_description(cls, description: str) -> 'MutationType':
        """Parse mutation type from description string."""
        if 'GATEKEEPER' in description:
            return cls.GATEKEEPER
        elif 'BOUNDARY' in description:
            return cls.BOUNDARY
        elif 'FLANKING' in description:
            return cls.FLANKING
        elif 'Direct mutation' in description:
            return cls.DIRECT
        elif 'Insertion' in description:
            return cls.INSERTION
        elif 'CORE' in description:
            return cls.BETA_CORE
        elif "Rule '" in description or 'MERGED RULES' in description:
            return cls.BETA_CORE
        return cls.UNKNOWN

    def __str__(self) -> str:
        return self.name


# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logging(
    verbose: bool = False,
    log_file: Optional[str] = None
) -> None:
    """
    Configure logging for the mutagenesis pipeline.
    
    Args:
        verbose: If True, show DEBUG level messages
        log_file: Optional path to write detailed logs
    """
    log_level = logging.DEBUG if verbose else logging.INFO

    formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%H:%M:%S'
    )

    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    # Clear existing handlers to avoid duplicates
    root_logger.handlers.clear()
    root_logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(funcName)s | %(message)s'
        )
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class PositionContext:
    """
    Structural context of a sequence position relative to aggregation-prone regions.
    
    Terminology:
    - APR: Aggregation-Prone Region (identified by clustering rules)
    - Gatekeeper: Position at APR boundary that can disrupt β-extension
    - Flanking: Beyond gatekeeper zone but still near APR
    - Core: Within APR cluster
    """
    position: int
    context_type: MutationType
    distance_to_nearest_apr: int  # 0 if within APR
    nearest_apr_boundary: Optional[int] = None
    apr_id: Optional[int] = None


@dataclass
class MutationInfo:
    """
    Structured representation of a mutation with biological context.
    
    Replaces ad-hoc dictionary construction with typed fields.
    """
    position: int
    original_aa: str
    new_aa: str
    mutation_type: MutationType
    aggregation_score: int
    region: Optional[str] = None
    rule_name: Optional[str] = None
    distance_to_apr_boundary: Optional[int] = None
    sequence: str = ""
    description: str = ""

    @property
    def is_gatekeeper(self) -> bool:
        return self.mutation_type == MutationType.GATEKEEPER

    @property
    def is_core(self) -> bool:
        return self.mutation_type == MutationType.BETA_CORE

    def to_tuple(self) -> Tuple[str, str, int]:
        """Convert to legacy (description, sequence, score) format."""
        return (self.description, self.sequence, self.aggregation_score)


@dataclass
class MultiMutationResult:
    """
    Structured result for multi-point mutations.
    
    Provides typed access to mutation combination properties and
    efficient classification methods.
    """
    description: str
    sequence: str
    aggregation_score: int
    positions: List[int]
    regions: List[str]
    level: int
    mutation_types: Tuple[MutationType, ...]
    type_composition: Tuple[str, ...]
    original_aas: Tuple[str, ...] = ()
    new_aas: Tuple[str, ...] = ()
    
    @property
    def is_all_gatekeeper(self) -> bool:
        """Check if all component mutations are gatekeepers."""
        return all(mt == MutationType.GATEKEEPER for mt in self.mutation_types)
    
    @property
    def is_all_core(self) -> bool:
        """Check if all component mutations are in β-core."""
        return all(mt == MutationType.BETA_CORE for mt in self.mutation_types)
    
    @property
    def is_mixed(self) -> bool:
        """Check if mutations span different types."""
        return len(set(self.mutation_types)) > 1
    
    @property
    def is_single_region(self) -> bool:
        """Check if all mutations are in the same region."""
        return len(set(self.regions)) <= 1
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'description': self.description,
            'sequence': self.sequence,
            'agg_score': self.aggregation_score,
            'positions': self.positions,
            'regions': self.regions,
            'level': self.level,
            'mutation_types': self.mutation_types,
            'type_composition': self.type_composition,
            'is_all_gatekeeper': self.is_all_gatekeeper,
            'is_all_core': self.is_all_core,
            'is_mixed': self.is_mixed,
        }
    
    def to_tuple(self) -> Tuple[str, str, int]:
        """Convert to legacy (description, sequence, score) format."""
        return (self.description, self.sequence, self.aggregation_score)


@dataclass(frozen=True, slots=True)
class Cluster:
    """
    Immutable representation of an aggregation-prone cluster.
    
    Using frozen=True enables:
    - Hashability (can be used in sets, as dict keys)
    - Thread safety (no mutation possible)
    - Semantic clarity (clusters are identified, not modified)
    
    Using slots=True reduces memory footprint by ~40% per instance.
    """
    positions: Tuple[int, ...]
    residues: Tuple[str, ...]
    rule_name: str
    aggregation_score: int

    def __post_init__(self):
        """Validate cluster invariants after initialization."""
        if len(self.positions) != len(self.residues):
            raise ValueError(
                f"Position count ({len(self.positions)}) must match "
                f"residue count ({len(self.residues)})"
            )
        if not self.positions:
            raise ValueError("Cluster must contain at least one position")

        if any(p < 1 for p in self.positions):
            raise ValueError("Positions must be 1-indexed (≥1)")
        if list(self.positions) != sorted(self.positions):
            raise ValueError("Positions must be in sorted order")

    @property
    def size(self) -> int:
        """Number of residues in the cluster."""
        return len(self.positions)

    @property
    def span(self) -> int:
        """Sequence span from first to last position (inclusive)."""
        return max(self.positions) - min(self.positions) + 1

    @property
    def density(self) -> float:
        """Cluster density: occupied positions / total span."""
        return self.size / self.span if self.span > 0 else 0.0

    @property
    def motif(self) -> str:
        """The cluster as a contiguous sequence string."""
        return ''.join(self.residues)

    @property
    def mean_propensity(self) -> float:
        """Average intrinsic aggregation propensity of cluster residues."""
        propensities = [AGGREGATION_PROPENSITY.get(aa, 0.0) for aa in self.residues]
        return sum(propensities) / len(propensities) if propensities else 0.0

    @property
    def min_position(self) -> int:
        """N-terminal boundary of cluster."""
        return min(self.positions)

    @property
    def max_position(self) -> int:
        """C-terminal boundary of cluster."""
        return max(self.positions)

    def overlaps_with(self, other: 'Cluster', gap_tolerance: int = 0) -> bool:
        """Check if this cluster overlaps with another."""
        self_min, self_max = self.min_position, self.max_position
        other_min, other_max = min(other.positions), max(other.positions)

        return not (
            self_max + gap_tolerance < other_min or
            other_max + gap_tolerance < self_min
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for backward compatibility."""
        return {
            'positions': list(self.positions),
            'residues': list(self.residues),
            'size': self.size,
            'span': self.span,
            'rule_name': self.rule_name,
            'aggregation_score': self.aggregation_score,
            'density': self.density,
            'mean_propensity': self.mean_propensity,
        }


@dataclass(frozen=True, slots=True)
class HydrophobicAromaticCluster:
    """
    Specialized cluster for hydrophobic-aromatic interaction patterns.
    
    Stores additional information about pair interactions and nearby
    hydrophobic residues, including sequence distances for validation
    of van der Waals contact feasibility.
    """
    positions: Tuple[int, ...]
    residues: Tuple[str, ...]
    rule_name: str
    aggregation_score: int
    pair_count: int = 0
    condition: str = ""  # 'at_least_2_pairs' or '1_pair_plus_hydrophobic'
    hydrophobic_count: int = 0
    aromatic_count: int = 0
    pair_interactions: Tuple[str, ...] = ()  # e.g., ("V10-F11", "L13-Y14")
    nearby_hydrophobic_positions: Tuple[int, ...] = ()
    nearby_hydrophobic_distances: Tuple[int, ...] = ()  # Parallel tuple of distances

    @property
    def size(self) -> int:
        return len(self.positions)

    @property
    def span(self) -> int:
        return max(self.positions) - min(self.positions) + 1 if self.positions else 0

    @property
    def min_position(self) -> int:
        return min(self.positions) if self.positions else 0

    @property
    def max_position(self) -> int:
        return max(self.positions) if self.positions else 0

    def overlaps_with(self, other: 'Cluster', gap_tolerance: int = 0) -> bool:
        """Check overlap with another cluster."""
        if not self.positions:
            return False
        self_min, self_max = self.min_position, self.max_position
        other_min, other_max = min(other.positions), max(other.positions)
        return not (
            self_max + gap_tolerance < other_min or
            other_max + gap_tolerance < self_min
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary with distance information preserved."""
        # Build nearby_hydrophobics list with distance and residue info
        nearby_hydrophobics = []
        for i, pos in enumerate(self.nearby_hydrophobic_positions):
            distance = (
                self.nearby_hydrophobic_distances[i]
                if i < len(self.nearby_hydrophobic_distances)
                else 0
            )
            # Get residue if position is in our positions tuple
            residue = '?'
            if pos in self.positions:
                idx = self.positions.index(pos)
                residue = self.residues[idx]
            nearby_hydrophobics.append({
                'position': pos,
                'distance': distance,
                'residue': residue
            })

        return {
            'positions': list(self.positions),
            'residues': list(self.residues),
            'size': self.size,
            'span': self.span,
            'rule_name': self.rule_name,
            'aggregation_score': self.aggregation_score,
            'pair_count': self.pair_count,
            'condition': self.condition,
            'hydrophobic_count': self.hydrophobic_count,
            'aromatic_count': self.aromatic_count,
            'pairs': [{'interaction': p} for p in self.pair_interactions],
            'nearby_hydrophobics': nearby_hydrophobics,
        }


@dataclass(frozen=True, slots=True)
class MultiRuleCluster:
    """
    Represents a cluster matching multiple aggregation rules.
    
    These clusters represent regions with compounded aggregation risk
    where multiple physicochemical features converge.
    """
    positions: Tuple[int, ...]
    residues: Tuple[str, ...]
    rules: Tuple[str, ...]
    combined_aggregation_score: int

    def __post_init__(self):
        if len(self.positions) != len(self.residues):
            raise ValueError("Position and residue counts must match")
        if len(self.rules) < 2:
            raise ValueError("MultiRuleCluster requires at least 2 rules")

    @property
    def size(self) -> int:
        return len(self.positions)

    @property
    def span(self) -> int:
        return max(self.positions) - min(self.positions) + 1 if self.positions else 0

    @property
    def min_position(self) -> int:
        return min(self.positions) if self.positions else 0

    @property
    def max_position(self) -> int:
        return max(self.positions) if self.positions else 0

    @property
    def is_multi_rule(self) -> bool:
        """Always True for MultiRuleCluster; aids duck typing."""
        return True

    @property
    def rule_string(self) -> str:
        """Joined rule names for display."""
        return '+'.join(sorted(self.rules))

    @property
    def aggregation_score(self) -> int:
        """Alias for combined_aggregation_score for duck typing."""
        return self.combined_aggregation_score

    @property
    def rule_name(self) -> str:
        """Alias for rule_string for duck typing."""
        return self.rule_string

    def overlaps_with(self, other: 'Cluster', gap_tolerance: int = 0) -> bool:
        """Check overlap with another cluster."""
        if not self.positions:
            return False
        self_min, self_max = self.min_position, self.max_position
        other_min, other_max = min(other.positions), max(other.positions)
        return not (
            self_max + gap_tolerance < other_min or
            other_max + gap_tolerance < self_min
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            'positions': list(self.positions),
            'residues': list(self.residues),
            'size': self.size,
            'span': self.span,
            'rules': list(self.rules),
            'combined_aggregation_score': self.combined_aggregation_score,
            'is_multi_rule': True,
        }


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
    
    Args:
        position: 1-indexed sequence position
        clusters: List of cluster dictionaries from analyze_region()
        gatekeeper_distance: Max distance from APR boundary for gatekeeper
        flanking_distance: Max distance for flanking region
    
    Returns:
        PositionContext with classification and distance metrics
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

        # Check if position is within cluster
        if position in positions:
            return PositionContext(
                position=position,
                context_type=MutationType.BETA_CORE,
                distance_to_nearest_apr=0,
                nearest_apr_boundary=position,
                apr_id=i,
            )

        # Calculate distance to cluster boundaries
        if position < cluster_min:
            distance = cluster_min - position
            boundary = cluster_min
        elif position > cluster_max:
            distance = position - cluster_max
            boundary = cluster_max
        else:
            # Position is between cluster positions but not IN cluster
            distance = min(abs(position - p) for p in positions)
            boundary = min(positions, key=lambda p: abs(position - p))

        if distance < min_distance:
            min_distance = distance
            nearest_boundary = boundary
            containing_cluster_id = i

    # Classify based on distance
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
        apr_id=containing_cluster_id if min_distance <= flanking_distance else None,
    )


def classify_mutation_type(
    position: int,
    new_aa: str,
    cluster_positions: List[int],
    cluster_min: int,
    cluster_max: int,
    gatekeeper_aas: Set[str] = None
) -> MutationType:
    """
    Determine mutation type based on position context AND amino acid properties.
    
    A true gatekeeper mutation requires BOTH:
    1. Position at APR boundary (within gatekeeper_distance)
    2. Mutation TO a gatekeeper amino acid (P, K, R, D, E)
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


# =============================================================================
# PROTOCOL AND EVALUATOR CLASSES
# =============================================================================

@runtime_checkable
class ClusterEvaluator(Protocol):
    """Protocol defining the interface for aggregation rule evaluators."""

    @property
    def name(self) -> str:
        ...

    @property
    def description(self) -> str:
        ...

    @property
    def aggregation_score(self) -> int:
        ...

    @property
    def residues(self) -> FrozenSet[str]:
        ...

    def find_clusters(
        self,
        sequence: str,
        start: int,
        stop: int
    ) -> List[Cluster]:
        ...


class BaseClusterEvaluator:
    """Base implementation providing common clustering logic."""

    NAME: str = "base"
    DESCRIPTION: str = "Base evaluator"
    RESIDUES: frozenset = frozenset()
    MIN_CLUSTER_SIZE: int = 2
    MAX_GAP: int = 3
    AGGREGATION_SCORE: int = 1

    @property
    def name(self) -> str:
        return self.NAME

    @property
    def description(self) -> str:
        return self.DESCRIPTION

    @property
    def aggregation_score(self) -> int:
        return self.AGGREGATION_SCORE

    @property
    def residues(self) -> frozenset:
        return self.RESIDUES

    def find_clusters(
        self,
        sequence: str,
        start: int,
        stop: int
    ) -> List[Cluster]:
        """Standard clustering implementation."""
        region_seq = sequence[start - 1:stop]

        matching_positions = [
            start + i
            for i, aa in enumerate(region_seq)
            if aa in self.RESIDUES
        ]

        raw_clusters = self._cluster_positions(matching_positions)

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
        """Group positions into clusters based on MAX_GAP threshold."""
        if not positions:
            return []

        clusters = []
        current = [positions[0]]

        for pos in positions[1:]:
            if pos - current[-1] <= self.MAX_GAP:
                current.append(pos)
            else:
                if len(current) >= 2:
                    clusters.append(current)
                current = [pos]

        if len(current) >= 2:
            clusters.append(current)

        return clusters


class HydrophobicAliphaticEvaluator(BaseClusterEvaluator):
    """Detects clusters of hydrophobic aliphatic residues (V, I, L, A, M)."""
    NAME = "hydrophobic_aliphatic"
    DESCRIPTION = "V, I, L, A, M clustered (≥3 residues within 4 positions)"
    RESIDUES = frozenset('VILAM')
    MIN_CLUSTER_SIZE = 3
    MAX_GAP = 4
    AGGREGATION_SCORE = 3


class AromaticEvaluator(BaseClusterEvaluator):
    """Detects clusters of aromatic residues (F, Y, W)."""
    NAME = "aromatic"
    DESCRIPTION = "F, Y, W clustered (≥2 residues within 3 positions)"
    RESIDUES = frozenset('FYW')
    MIN_CLUSTER_SIZE = 2
    MAX_GAP = 3
    AGGREGATION_SCORE = 2


class AmideEvaluator(BaseClusterEvaluator):
    """Detects clusters of amide-containing residues (Q, N)."""
    NAME = "amide"
    DESCRIPTION = "Q, N clustered (≥2 residues within 3 positions)"
    RESIDUES = frozenset('QN')
    MIN_CLUSTER_SIZE = 2
    MAX_GAP = 3
    AGGREGATION_SCORE = 1


class HydrophobicAromaticEvaluator(BaseClusterEvaluator):
    """
    Detects hydrophobic-aromatic adjacency patterns.
    
    Special logic:
    - Condition A: ≥2 hydrophobic-aromatic adjacent pairs
    - Condition B: 1 pair + additional hydrophobic within 3 positions
    """
    NAME = "hydrophobic_and_aromatic"
    DESCRIPTION = "Hydrophobic (V,I,L,A,M) adjacent to aromatic (F,Y,W)"
    RESIDUES = frozenset('VILAMFYW')
    MIN_CLUSTER_SIZE = 2
    MAX_GAP = 1
    AGGREGATION_SCORE = 2

    HYDROPHOBIC = frozenset('VILAM')
    AROMATIC = frozenset('FYW')

    def find_clusters(
        self,
        sequence: str,
        start: int,
        stop: int
    ) -> List[Union[Cluster, HydrophobicAromaticCluster]]:
        """Override base implementation with special pair-detection logic."""
        region_seq = sequence[start - 1:stop]

        # Find adjacent hydrophobic-aromatic pairs
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

        # Collect all hydrophobic positions for condition B
        hydrophobic_positions = {
            start + i
            for i, aa in enumerate(region_seq)
            if aa in self.HYDROPHOBIC
        }

        clusters = []

        # Condition A: ≥2 pairs close together
        if len(pairs) >= 2:
            pair_positions = set()
            pair_interactions = []
            for pair in pairs:
                pair_positions.update(pair['positions'])
                pair_interactions.append(pair['interaction'])

            sorted_positions = tuple(sorted(pair_positions))
            cluster_residues = tuple(sequence[p - 1] for p in sorted_positions)

            clusters.append(HydrophobicAromaticCluster(
                positions=sorted_positions,
                residues=cluster_residues,
                rule_name=self.NAME,
                aggregation_score=self.AGGREGATION_SCORE,
                pair_count=len(pairs),
                condition='at_least_2_pairs',
                hydrophobic_count=sum(1 for r in cluster_residues if r in self.HYDROPHOBIC),
                aromatic_count=sum(1 for r in cluster_residues if r in self.AROMATIC),
                pair_interactions=tuple(pair_interactions),
                nearby_hydrophobic_positions=(),
                nearby_hydrophobic_distances=()
            ))

        # Condition B: 1 pair + nearby hydrophobic
        for pair in pairs:
            pair_position_set = set(pair['positions'])

            nearby = []
            nearby_distances = []
            for hp in hydrophobic_positions:
                if hp in pair_position_set:
                    continue
                min_distance = min(abs(hp - p) for p in pair['positions'])
                if min_distance <= 3:
                    nearby.append(hp)
                    nearby_distances.append(min_distance)

            if nearby:
                all_positions = tuple(sorted(pair_position_set | set(nearby)))
                cluster_residues = tuple(sequence[p - 1] for p in all_positions)

                clusters.append(HydrophobicAromaticCluster(
                    positions=all_positions,
                    residues=cluster_residues,
                    rule_name=self.NAME,
                    aggregation_score=self.AGGREGATION_SCORE,
                    pair_count=1,
                    condition='1_pair_plus_hydrophobic',
                    hydrophobic_count=sum(1 for r in cluster_residues if r in self.HYDROPHOBIC),
                    aromatic_count=sum(1 for r in cluster_residues if r in self.AROMATIC),
                    pair_interactions=(pair['interaction'],),
                    nearby_hydrophobic_positions=tuple(nearby),
                    nearby_hydrophobic_distances=tuple(nearby_distances)
                ))

        # Deduplicate by positions
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
    """Registry for managing cluster evaluation rules."""

    def __init__(self):
        self._evaluators: Dict[str, ClusterEvaluator] = {}

    def register(self, evaluator: ClusterEvaluator) -> None:
        """Register an evaluator under its name."""
        if evaluator.name in self._evaluators:
            logger.warning(f"Overwriting existing evaluator: {evaluator.name}")
        self._evaluators[evaluator.name] = evaluator
        logger.debug(f"Registered evaluator: {evaluator.name}")

    def get(self, name: str) -> ClusterEvaluator:
        """Retrieve evaluator by name."""
        if name not in self._evaluators:
            raise KeyError(f"Unknown rule: {name}. Available: {list(self._evaluators.keys())}")
        return self._evaluators[name]

    def list_rules(self) -> List[str]:
        """Get list of all registered rule names."""
        return list(self._evaluators.keys())

    def evaluate_region(
        self,
        sequence: str,
        start: int,
        stop: int,
        rules: Optional[List[str]] = None
    ) -> Dict[str, List[Cluster]]:
        """Evaluate all (or selected) rules on a region."""
        target_rules = rules or self.list_rules()

        results = {}
        for rule_name in target_rules:
            if rule_name not in self._evaluators:
                logger.warning(f"Skipping unknown rule: {rule_name}")
                continue

            evaluator = self._evaluators[rule_name]
            clusters = evaluator.find_clusters(sequence, start, stop)
            results[rule_name] = clusters

            logger.debug(
                f"Rule '{rule_name}' found {len(clusters)} clusters in {start}:{stop}"
            )

        return results


def create_default_registry() -> RuleRegistry:
    """Create registry with all standard aggregation rules."""
    registry = RuleRegistry()
    registry.register(HydrophobicAliphaticEvaluator())
    registry.register(AromaticEvaluator())
    registry.register(AmideEvaluator())
    registry.register(HydrophobicAromaticEvaluator())
    return registry


# Module-level default registry
DEFAULT_REGISTRY = create_default_registry()


def _generate_rules_dict(registry: RuleRegistry) -> Dict[str, Dict]:
    """Generate RULES dictionary from registry for backward compatibility."""
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


RULES = _generate_rules_dict(DEFAULT_REGISTRY)


# =============================================================================
# UNION-FIND DATA STRUCTURE
# =============================================================================

class UnionFind:
    """Disjoint Set Union (Union-Find) with path compression and union by rank."""

    __slots__ = ('parent', 'rank', '_size')

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n
        self._size = [1] * n

    def find(self, x: int) -> int:
        """Find root with path compression."""
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x: int, y: int) -> bool:
        """Merge sets with union by rank."""
        px, py = self.find(x), self.find(y)

        if px == py:
            return False

        if self.rank[px] < self.rank[py]:
            px, py = py, px

        self.parent[py] = px
        self._size[px] += self._size[py]

        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1

        return True

    def connected(self, x: int, y: int) -> bool:
        """Check if x and y are in the same set."""
        return self.find(x) == self.find(y)

    def set_size(self, x: int) -> int:
        """Get the size of the set containing x."""
        return self._size[self.find(x)]

    def get_groups(self) -> Dict[int, List[int]]:
        """Get all disjoint sets as a dictionary."""
        groups: Dict[int, List[int]] = {}
        for i in range(len(self.parent)):
            root = self.find(i)
            if root not in groups:
                groups[root] = []
            groups[root].append(i)
        return groups


# =============================================================================
# CLUSTER OPERATIONS
# =============================================================================

def create_cluster(
    positions: List[int],
    sequence: str,
    rule_name: str,
    aggregation_score: int
) -> Cluster:
    """Factory function to create Cluster from mutable inputs."""
    sorted_positions = tuple(sorted(positions))
    residues = tuple(sequence[p - 1] for p in sorted_positions)

    return Cluster(
        positions=sorted_positions,
        residues=residues,
        rule_name=rule_name,
        aggregation_score=aggregation_score
    )


def merge_overlapping_clusters_unionfind(
    clusters: List[Cluster],
    sequence: str,
    gap_tolerance: int = MAX_GAP_FOR_MERGING
) -> List[Union[Cluster, MultiRuleCluster]]:
    """Merge overlapping or proximal clusters using Union-Find."""
    if not clusters:
        return []

    n = len(clusters)
    uf = UnionFind(n)

    # Phase 1: Identify all overlapping pairs and union them
    for i in range(n):
        for j in range(i + 1, n):
            if clusters[i].overlaps_with(clusters[j], gap_tolerance):
                uf.union(i, j)

    # Phase 2: Group clusters by their root
    groups = uf.get_groups()

    # Phase 3: Merge each group into a single cluster
    merged_clusters = []

    for root, indices in groups.items():
        group_clusters = [clusters[i] for i in indices]

        if len(group_clusters) == 1:
            merged_clusters.append(group_clusters[0])
        else:
            merged = _merge_cluster_group_to_object(group_clusters, sequence)
            merged_clusters.append(merged)

    return merged_clusters


def _merge_cluster_group_to_object(
    clusters: List[Cluster],
    sequence: str
) -> Union[Cluster, MultiRuleCluster]:
    """Merge a group of overlapping clusters into a single object."""
    all_positions: Set[int] = set()
    all_rules: Set[str] = set()
    max_score = 0
    total_score = 0

    for cluster in clusters:
        all_positions.update(cluster.positions)
        all_rules.add(cluster.rule_name)
        max_score = max(max_score, cluster.aggregation_score)
        total_score += cluster.aggregation_score

    sorted_positions = tuple(sorted(all_positions))
    residues = tuple(sequence[p - 1] for p in sorted_positions)

    if len(all_rules) > 1:
        combined_score = total_score + len(all_rules) - 1
        return MultiRuleCluster(
            positions=sorted_positions,
            residues=residues,
            rules=tuple(sorted(all_rules)),
            combined_aggregation_score=combined_score
        )
    else:
        rule_name = next(iter(all_rules))
        return Cluster(
            positions=sorted_positions,
            residues=residues,
            rule_name=rule_name,
            aggregation_score=max_score
        )


def find_clusters(positions: List[int], max_gap: int = 3) -> List[List[int]]:
    """Find clusters of positions where positions are within max_gap of each other."""
    if not positions:
        return []

    clusters = []
    current_cluster = [positions[0]]

    for i in range(1, len(positions)):
        if positions[i] - current_cluster[-1] <= max_gap:
            current_cluster.append(positions[i])
        else:
            if len(current_cluster) >= 2:
                clusters.append(current_cluster)
            current_cluster = [positions[i]]

    if len(current_cluster) >= 2:
        clusters.append(current_cluster)

    return clusters


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def analyze_region(
    sequence: str,
    start: int,
    stop: int,
    registry: RuleRegistry = None,
    selected_rules: Optional[List[str]] = None
) -> Dict[str, Any]:
    """Analyze a protein region for aggregation-prone clusters."""
    if registry is None:
        registry = DEFAULT_REGISTRY

    region_seq = sequence[start - 1:stop]

    rule_results = registry.evaluate_region(sequence, start, stop, selected_rules)

    all_clusters: List[Cluster] = []
    for rule_name, clusters in rule_results.items():
        all_clusters.extend(clusters)

    merged_clusters = merge_overlapping_clusters_unionfind(
        all_clusters, sequence, gap_tolerance=MAX_GAP_FOR_MERGING
    )

    multi_rule_clusters = [
        c for c in merged_clusters
        if isinstance(c, MultiRuleCluster)
    ]

    hotspot_positions = set()
    for cluster in merged_clusters:
        hotspot_positions.update(cluster.positions)

    results = {
        'region': (start, stop),
        'sequence': region_seq,
        'length': len(region_seq),
        'rules': {},
        'merged_clusters': [c.to_dict() for c in merged_clusters],
        'multi_rule_clusters': [c.to_dict() for c in multi_rule_clusters],
        'aggregation_hotspots': sorted(hotspot_positions),
    }

    for rule_name, clusters in rule_results.items():
        evaluator = registry.get(rule_name)

        rule_entry = {
            'description': evaluator.description,
            'qualifying_clusters': [c.to_dict() for c in clusters],
            'condition_met': len(clusters) > 0,
            'aggregation_score': evaluator.aggregation_score,
            'matching_positions': [],
            'matching_residues': [],
            'clusters': [],
        }

        for cluster in clusters:
            rule_entry['matching_positions'].extend(cluster.positions)
            rule_entry['matching_residues'].extend(cluster.residues)
            rule_entry['clusters'].append(list(cluster.positions))

        rule_entry['matching_positions'] = list(dict.fromkeys(rule_entry['matching_positions']))
        rule_entry['matching_residues'] = list(dict.fromkeys(rule_entry['matching_residues']))

        if rule_name == 'hydrophobic_and_aromatic' and clusters:
            rule_entry['special_clusters'] = [c.to_dict() for c in clusters]

        results['rules'][rule_name] = rule_entry

    return results


# =============================================================================
# MUTATION EXTRACTION AND PARSING
# =============================================================================

def extract_mutation_info(
    description: str,
    sequence: str = "",
    regions: Optional[List[str]] = None
) -> MutationInfo:
    """
    Extract structured mutation information from description string.
    
    Uses explicit patterns to avoid matching aggregation scores as positions.
    """
    # Pattern for standard mutations: OriginalAA + Position + NewAA at start
    mutation_pattern = re.compile(r'^([A-Z])(\d+)([A-Z])\s*\|')

    # Pattern for insertions
    insertion_pattern = re.compile(
        r'Insertion:\s*([A-Z])\s+inserted\s+before\s+position\s+(\d+)'
    )

    # Pattern for aggregation score
    score_pattern = re.compile(r'agg_score=(\d+)')

    # Pattern for rule name
    rule_pattern = re.compile(r"Rule\s+'([^']+)'")
    merged_pattern = re.compile(r'MERGED RULES\s+([\w+]+)')

    # Initialize defaults
    position = None
    original_aa = None
    new_aa = None
    agg_score = 0
    rule_name = None
    mutation_type = MutationType.UNKNOWN

    # Try standard mutation format
    match = mutation_pattern.match(description)
    if match:
        original_aa = match.group(1)
        position = int(match.group(2))
        new_aa = match.group(3)
    else:
        # Try insertion format
        ins_match = insertion_pattern.search(description)
        if ins_match:
            new_aa = ins_match.group(1)
            position = int(ins_match.group(2))
            original_aa = '-'
            mutation_type = MutationType.INSERTION

    # Extract score
    score_match = score_pattern.search(description)
    if score_match:
        agg_score = int(score_match.group(1))

    # Extract rule
    rule_match = rule_pattern.search(description)
    if rule_match:
        rule_name = rule_match.group(1)
    else:
        merged_match = merged_pattern.search(description)
        if merged_match:
            rule_name = merged_match.group(1)

    # Determine mutation type from description
    if mutation_type == MutationType.UNKNOWN:
        mutation_type = MutationType.from_description(description)

    # Determine region
    region = None
    if regions and position:
        for region_str in regions:
            try:
                start, stop = map(int, region_str.split(':'))
                if start <= position <= stop:
                    region = region_str
                    break
            except ValueError:
                continue

    return MutationInfo(
        position=position or 0,
        original_aa=original_aa or '?',
        new_aa=new_aa or '?',
        mutation_type=mutation_type,
        aggregation_score=agg_score,
        region=region,
        rule_name=rule_name,
        description=description,
        sequence=sequence,
    )


def get_region_for_position(
    position: int,
    regions: Optional[List[str]] = None
) -> Optional[str]:
    """Get the region(s) a position belongs to."""
    if not regions or position is None:
        return None

    matching_regions = []
    for region_str in regions:
        try:
            start, stop = map(int, region_str.split(':'))
            if start <= position <= stop:
                matching_regions.append(region_str)
        except ValueError:
            continue

    return ','.join(matching_regions) if matching_regions else None


# =============================================================================
# MUTATION APPLICATION
# =============================================================================

def create_mutated_sequence(sequence: str, position: int, new_aa: str) -> str:
    """Create a mutated sequence by replacing one amino acid."""
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} out of range (1-{len(sequence)})")

    if len(sequence) < 1000:
        seq_list = list(sequence)
        seq_list[position - 1] = new_aa
        return ''.join(seq_list)
    else:
        return sequence[:position-1] + new_aa + sequence[position:]


def get_mutations_for_position(
    pos: int,
    cluster_positions: List[int],
    region_start: int,
    region_end: int,
    mutations: List[str],
    gatekeeping_aas: List[str]
) -> List[str]:
    """Get list of mutations to apply based on position type."""
    if is_edge_position(pos, cluster_positions, region_start, region_end):
        return list(set(mutations + gatekeeping_aas))
    return mutations


def is_edge_position(
    pos: int,
    cluster_positions: List[int],
    region_start: int,
    region_end: int
) -> bool:
    """Check if a position is at the edge of a cluster."""
    if not cluster_positions:
        return False

    min_pos = min(cluster_positions)
    max_pos = max(cluster_positions)

    if pos == min_pos or pos == max_pos:
        return True

    if pos == region_start or pos == region_end:
        return True

    if pos == region_start + 1 or pos == region_end - 1:
        return True

    return False


def apply_rule_mutations(
    sequence: str,
    region_analysis: Dict,
    mutations: List[str],
    selected_rules: List[str] = None,
    gatekeeping_aas: List[str] = None,
    gatekeeper_distance: int = GATEKEEPER_DISTANCE
) -> List[Tuple[str, str, int, MutationType]]:
    """
    Apply mutations based on rules with biologically accurate classification.
    
    Returns tuples of (description, sequence, agg_score, mutation_type).
    """
    results = []

    if gatekeeping_aas is None:
        gatekeeping_aas = GATEKEEPING_AAS

    gatekeeper_aa_set = set(aa.upper() for aa in gatekeeping_aas) | CANONICAL_GATEKEEPER_AAS

    # Filter rules
    if selected_rules:
        rules_to_apply = {
            k: v for k, v in region_analysis['rules'].items()
            if k in selected_rules
        }
    else:
        rules_to_apply = region_analysis['rules']

    region_start, region_end = region_analysis['region']

    # Collect all clusters from all rules
    cluster_objects = []
    for rule_name, rule_data in rules_to_apply.items():
        if not rule_data['condition_met']:
            continue

        for cluster_dict in rule_data['qualifying_clusters']:
            cluster_objects.append(create_cluster(
                positions=cluster_dict['positions'],
                sequence=sequence,
                rule_name=cluster_dict.get('rule_name', rule_name),
                aggregation_score=cluster_dict.get('aggregation_score', rule_data['aggregation_score'])
            ))

    # Merge overlapping clusters
    merged_clusters = merge_overlapping_clusters_unionfind(
        cluster_objects, sequence, gap_tolerance=MAX_GAP_FOR_MERGING
    )

    # Track mutated positions
    mutated_positions: Set[int] = set()

    # Apply mutations to each merged cluster
    for cluster in merged_clusters:
        positions_to_mutate = list(cluster.positions)

        # Get cluster boundaries
        cluster_min = min(positions_to_mutate)
        cluster_max = max(positions_to_mutate)

        # Get aggregation score and rule description
        if isinstance(cluster, MultiRuleCluster):
            agg_score = cluster.combined_aggregation_score
            rule_desc = f"MERGED RULES {'+'.join(cluster.rules)}"
        else:
            agg_score = cluster.aggregation_score
            rule_desc = f"Rule '{cluster.rule_name}'"

        for pos in positions_to_mutate:
            if pos in mutated_positions:
                continue

            mutated_positions.add(pos)
            original_aa = sequence[pos - 1]

            # Determine if this is a boundary position
            is_boundary = (pos == cluster_min or pos == cluster_max)

            # Get applicable mutations
            if is_boundary:
                all_mutations = list(set(mutations + gatekeeping_aas))
            else:
                all_mutations = mutations

            for new_aa in all_mutations:
                if new_aa == original_aa:
                    continue

                mutated_seq = create_mutated_sequence(sequence, pos, new_aa)

                # Classify mutation type
                mutation_type = classify_mutation_type(
                    position=pos,
                    new_aa=new_aa,
                    cluster_positions=positions_to_mutate,
                    cluster_min=cluster_min,
                    cluster_max=cluster_max,
                    gatekeeper_aas=gatekeeper_aa_set
                )

                # Build description with type label
                type_label = mutation_type.name
                description = (
                    f"{original_aa}{pos}{new_aa} | {rule_desc} | "
                    f"{type_label} (agg_score={agg_score})"
                )

                results.append((description, mutated_seq, agg_score, mutation_type))

    return results


def mutate_sequence(
    sequence: str,
    positions: List[int],
    mutations: List[str],
    regions: List[str] = None,
    selected_rules: List[str] = None,
    insertion_positions: List[int] = None,
    insertion_aas: List[str] = None,
    gatekeeping_aas: List[str] = None,
    verbose: bool = False
) -> Tuple[List[Tuple[str, str]], List[Dict]]:
    """Generate mutated sequences with point mutations, rule-based mutations, and insertions."""
    all_results = []
    region_analyses = []
    seq_len = len(sequence)

    if gatekeeping_aas is None:
        gatekeeping_aas = GATEKEEPING_AAS

    # Validate positions
    for pos in positions:
        if pos < 1 or pos > seq_len:
            raise ValueError(f"Position {pos} is out of range (sequence length: {seq_len})")

    # Analyze regions for rule-based mutations
    if regions:
        for region_str in regions:
            start, stop = parse_region(region_str, seq_len)
            analysis = analyze_region(sequence, start, stop)
            region_analyses.append(analysis)

            logger.debug(f"Analyzing region {start}:{stop}")
            logger.debug(f"Sequence: {analysis['sequence']}")
            logger.debug(f"Length: {analysis['length']} residues")

            # Verbose logging for rules
            for rule_name, rule_data in analysis['rules'].items():
                if rule_name == 'hydrophobic_and_aromatic' and rule_data['condition_met']:
                    logger.debug(f"\n{rule_name.upper()}:")
                    logger.debug(f"  Description: {rule_data['description']}")
                    logger.debug(f"  ✓ Rule triggered!")

                    for i, cluster in enumerate(rule_data.get('special_clusters', []), 1):
                        logger.debug(f"\n  Cluster {i} ({cluster['condition']}):")
                        logger.debug(f"    Positions: {cluster['positions']}")
                        logger.debug(f"    Residues: {''.join(cluster['residues'])}")
                        logger.debug(f"    Hydrophobic count: {cluster['hydrophobic_count']}")
                        logger.debug(f"    Aromatic count: {cluster['aromatic_count']}")

                        if cluster['condition'] == 'at_least_2_pairs':
                            logger.debug(f"    Pairs found: {cluster['pair_count']}")
                            for j, pair in enumerate(cluster.get('pairs', []), 1):
                                logger.debug(f"      Pair {j}: {pair['interaction']}")
                        else:
                            if cluster.get('pairs'):
                                pair = cluster['pairs'][0]
                                logger.debug(f"    Pair: {pair['interaction']}")
                            if cluster.get('nearby_hydrophobics'):
                                logger.debug(f"    Nearby hydrophobic residues:")
                                for h in cluster['nearby_hydrophobics']:
                                    pos = h.get('position', 'unknown')
                                    dist = h.get('distance', '?')
                                    res = h.get('residue', '?')
                                    logger.debug(f"      {res}{pos} (distance: {dist})")

                        logger.debug(f"    Size: {cluster['size']}, Span: {cluster['span']}aa")

                elif rule_data['matching_positions']:
                    logger.debug(f"\n{rule_name.upper()}:")
                    logger.debug(f"  Description: {rule_data['description']}")
                    residues_str = ''.join(rule_data['matching_residues'])
                    logger.debug(f"  Matching residues: {residues_str}")
                    logger.debug(f"  Positions: {rule_data['matching_positions']}")

                    if rule_data['clusters']:
                        logger.debug(f"  Found {len(rule_data['clusters'])} cluster(s):")
                        for i, cluster in enumerate(rule_data['clusters'], 1):
                            cluster_residues = [sequence[pos-1] for pos in cluster]
                            span = max(cluster) - min(cluster) + 1
                            logger.debug(
                                f"    Cluster {i}: positions {cluster}, "
                                f"residues {''.join(cluster_residues)}, "
                                f"size={len(cluster)}, span={span}aa"
                            )

                    if rule_data['qualifying_clusters']:
                        logger.debug(f"  ✓ QUALIFYING CLUSTERS:")
                        for i, cluster in enumerate(rule_data['qualifying_clusters'], 1):
                            logger.debug(
                                f"    Cluster {i}: positions {cluster['positions']}, "
                                f"residues {''.join(cluster['residues'])}, "
                                f"size={cluster['size']}, span={cluster['span']}aa, "
                                f"agg_score={cluster['aggregation_score']}"
                            )
                    else:
                        logger.debug(f"  ✗ No qualifying clusters")

            if analysis['multi_rule_clusters']:
                logger.debug(f"\n{'*'*60}")
                logger.debug(f"MULTI-RULE CLUSTERS (HIGH AGGREGATION RISK):")
                logger.debug(f"{'*'*60}")
                for i, cluster in enumerate(analysis['multi_rule_clusters'], 1):
                    logger.debug(f"\n  Multi-Rule Cluster {i}:")
                    logger.debug(f"    Rules: {', '.join(cluster['rules'])}")
                    logger.debug(f"    Positions: {cluster['positions']}")
                    logger.debug(f"    Residues: {''.join(cluster['residues'])}")
                    logger.debug(f"    Combined Aggregation Score: {cluster['combined_aggregation_score']}/8")

    # Apply rule-based mutations
    for analysis in region_analyses:
        rule_mutations = apply_rule_mutations(
            sequence, analysis, mutations, selected_rules, gatekeeping_aas
        )
        # Convert to format with mutation_type preserved in description
        for desc, seq, score, mut_type in rule_mutations:
            all_results.append((desc, seq, score, mut_type))

    # Generate direct point mutations
    for pos in positions:
        original_aa = sequence[pos-1]

        for new_aa in mutations:
            if new_aa == original_aa:
                continue

            mutated_seq = create_mutated_sequence(sequence, pos, new_aa)
            description = f"{original_aa}{pos}{new_aa} | Direct mutation | DIRECT (agg_score=0)"
            all_results.append((description, mutated_seq, 0, MutationType.DIRECT))

    # Generate insertions
    if insertion_positions and insertion_aas:
        if len(insertion_positions) != len(insertion_aas):
            raise ValueError("insertion_positions and insertion_aas must have the same length")

        for ins_pos, ins_aa in zip(insertion_positions, insertion_aas):
            if ins_pos < 1 or ins_pos > seq_len + 1:
                raise ValueError(f"Insertion position {ins_pos} is out of range (1-{seq_len+1})")

            mutated_seq = sequence[:ins_pos-1] + ins_aa + sequence[ins_pos-1:]
            description = f"Insertion: {ins_aa} inserted before position {ins_pos} | INSERTION (agg_score=0)"
            all_results.append((description, mutated_seq, 0, MutationType.INSERTION))

    # Sort by aggregation score (descending)
    all_results.sort(key=lambda x: x[2], reverse=True)

    # Convert to format expected by downstream functions
    sorted_results = [(desc, seq) for desc, seq, _, _ in all_results]

    return sorted_results, region_analyses


# =============================================================================
# MULTI-POINT MUTATION GENERATION
# =============================================================================

def _parse_mutation_infos(
    mutations: List[Tuple[str, str]],
    sequence: str,
    regions: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """Parse mutation descriptions into structured info dictionaries."""
    mutation_infos = []

    for desc, seq in mutations:
        info = extract_mutation_info(desc, seq, regions)
        if not info.position or not seq:
            continue

        # Only point mutations can be combined
        match = re.search(r'([A-Z])(\d+)([A-Z])', desc)
        if not match:
            continue

        mutation_infos.append({
            'position': info.position,
            'mutation_type': info.mutation_type,
            'is_gatekeeper': info.is_gatekeeper,
            'region': info.region,
            'agg_score': info.aggregation_score,
            'description': desc,
            'sequence': seq,
            'original_aa': sequence[info.position - 1] if info.position <= len(sequence) else '?',
            'new_aa': match.group(3)
        })

    return mutation_infos


def _filter_by_position(
    mutation_infos: List[Dict[str, Any]],
    top_variants_per_position: int = 3
) -> Dict[int, List[Dict[str, Any]]]:
    """Group and filter mutations by position."""
    by_position: Dict[int, List[Dict[str, Any]]] = {}

    for mi in mutation_infos:
        by_position.setdefault(mi['position'], []).append(mi)

    # Deduplicate per-position variants by new amino acid
    filtered: Dict[int, List[Dict[str, Any]]] = {}

    for pos, variants in by_position.items():
        best_by_newaa: Dict[str, Dict[str, Any]] = {}

        for v in variants:
            new_aa = v.get('new_aa')
            if not new_aa:
                continue
            existing = best_by_newaa.get(new_aa)
            if existing is None or v.get('agg_score', 0) > existing.get('agg_score', 0):
                best_by_newaa[new_aa] = v

        kept = list(best_by_newaa.values())
        kept.sort(key=lambda x: x.get('agg_score', 0), reverse=True)

        if top_variants_per_position and top_variants_per_position > 0:
            kept = kept[:top_variants_per_position]

        filtered[pos] = kept

    return filtered


def _create_multi_mutation_result(
    variant_combo: Tuple[Dict[str, Any], ...],
    sequence: str,
    level: int,
    regions: Optional[List[str]]
) -> MultiMutationResult:
    """Create a MultiMutationResult from a variant combination."""
    positions = [v['position'] for v in variant_combo]
    new_aas = tuple(v['new_aa'] for v in variant_combo)
    original_aas = tuple(v['original_aa'] for v in variant_combo)
    mutation_types = tuple(v.get('mutation_type', MutationType.UNKNOWN) for v in variant_combo)

    # Apply mutations in reverse order
    mutations_to_apply = sorted(zip(positions, new_aas), reverse=True, key=lambda x: x[0])
    combined_seq = sequence
    for pos, aa in mutations_to_apply:
        combined_seq = create_mutated_sequence(combined_seq, pos, aa)

    combined_agg_score = sum(v.get('agg_score', 0) for v in variant_combo)

    mutations_str = ' + '.join(
        f"{v['original_aa']}{v['position']}{v['new_aa']}"
        for v in sorted(variant_combo, key=lambda x: x['position'])
    )

    regions_involved = sorted({v['region'] for v in variant_combo if v.get('region')})
    regions_str = ', '.join(regions_involved) if regions_involved else "Direct"

    # Build type composition string
    type_composition = tuple(
        mt.name if isinstance(mt, MutationType) else str(mt)
        for mt in mutation_types
    )

    # Create result object
    result = MultiMutationResult(
        description="",  # Will be set below
        sequence=combined_seq,
        aggregation_score=combined_agg_score,
        positions=positions,
        regions=regions_involved,
        level=level,
        mutation_types=mutation_types,
        type_composition=type_composition,
        original_aas=original_aas,
        new_aas=new_aas,
    )

    # Build description based on classification
    description = f"{mutations_str} | {regions_str} (level={level}, agg_score={combined_agg_score})"
    if result.is_all_gatekeeper:
        description += " | ALL_GATEKEEPER"
    elif result.is_all_core:
        description += " | ALL_CORE"
    elif result.is_mixed:
        description += " | MIXED"

    # Update description (dataclass is not frozen so we can do this)
    object.__setattr__(result, 'description', description)

    return result


def _process_position_chunk(
    position_combos: List[Tuple[int, ...]],
    filtered_by_position: Dict[int, List[Dict]],
    sequence: str,
    level: int,
    regions: Optional[List[str]],
    max_per_chunk: int
) -> List[Dict[str, Any]]:
    """Process a chunk of position combinations (worker function)."""
    results = []
    count = 0

    for pos_combo in position_combos:
        if count >= max_per_chunk:
            break

        variant_lists = [filtered_by_position[p] for p in pos_combo]

        for variant_combo in product(*variant_lists):
            if count >= max_per_chunk:
                break

            result = _create_multi_mutation_result(variant_combo, sequence, level, regions)
            results.append(result.to_dict())
            count += 1

    return results


def generate_multi_point_mutations(
    mutations: List[Tuple[str, str]],
    sequence: str,
    multi_mutation_levels: Optional[List[int]] = None,
    regions: Optional[List[str]] = None,
    max_combinations: int = 10000,
    top_variants_per_position: int = 3,
    n_workers: Optional[int] = None
) -> Dict[int, List[Dict[str, Any]]]:
    """
    Generate multi-point mutations from single mutations.
    
    Supports parallel processing for large combination spaces.
    
    Args:
        mutations: List of (description, sequence) tuples for single mutations
        sequence: Original protein sequence
        multi_mutation_levels: Levels to generate (e.g., [2, 3])
        regions: List of region strings for categorization
        max_combinations: Maximum combinations per level
        top_variants_per_position: Limit variants per position
        n_workers: Number of parallel workers (None = auto)
    
    Returns:
        Dict mapping level to list of mutation result dicts
    """
    # Guard clause for empty/None input
    if not multi_mutation_levels:
        logger.warning("No multi-mutation levels specified; returning empty results")
        return {}

    # Validate levels
    valid_levels = []
    for level in multi_mutation_levels:
        if not isinstance(level, int):
            logger.warning(f"Skipping non-integer level: {level}")
            continue
        if level < 2:
            logger.warning(f"Skipping invalid level {level} (must be ≥2)")
            continue
        if level > MAX_MUTATION_LEVEL:
            logger.warning(f"Skipping level {level} (exceeds maximum {MAX_MUTATION_LEVEL})")
            continue
        valid_levels.append(level)

    if not valid_levels:
        logger.warning("No valid multi-mutation levels after filtering")
        return {}

    # Parse mutation info
    mutation_infos = _parse_mutation_infos(mutations, sequence, regions)
    filtered_by_position = _filter_by_position(mutation_infos, top_variants_per_position)
    viable_positions = sorted([p for p, vs in filtered_by_position.items() if vs])

    if not viable_positions:
        logger.warning("No viable positions for multi-mutations")
        return {}

    # Determine worker count
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    multi_results = {}

    for level in valid_levels:
        if level > len(viable_positions):
            logger.warning(f"Skipping level {level} (only {len(viable_positions)} positions available)")
            continue

        # Generate position combinations
        position_combos = list(combinations(viable_positions, level))

        if not position_combos:
            multi_results[level] = []
            continue

        logger.info(f"Level {level}: Processing {len(position_combos)} position combinations")

        # Decide whether to parallelize
        if len(position_combos) > 1000 and n_workers > 1:
            # Parallel processing
            n_chunks = min(n_workers, len(position_combos))
            chunk_size = len(position_combos) // n_chunks + 1
            chunks = [
                position_combos[i:i + chunk_size]
                for i in range(0, len(position_combos), chunk_size)
            ]

            all_results = []
            try:
                with ProcessPoolExecutor(max_workers=n_workers) as executor:
                    futures = {
                        executor.submit(
                            _process_position_chunk,
                            chunk,
                            filtered_by_position,
                            sequence,
                            level,
                            regions,
                            max_combinations // n_chunks
                        ): i
                        for i, chunk in enumerate(chunks)
                    }

                    for future in as_completed(futures):
                        try:
                            chunk_results = future.result()
                            all_results.extend(chunk_results)
                        except Exception as e:
                            logger.error(f"Chunk processing failed: {e}")

            except Exception as e:
                logger.warning(f"Parallel processing failed, falling back to sequential: {e}")
                all_results = _process_position_chunk(
                    position_combos,
                    filtered_by_position,
                    sequence,
                    level,
                    regions,
                    max_combinations
                )
        else:
            # Sequential processing
            all_results = _process_position_chunk(
                position_combos,
                filtered_by_position,
                sequence,
                level,
                regions,
                max_combinations
            )

        # Sort and limit
        all_results.sort(key=lambda x: x.get('agg_score', 0), reverse=True)
        multi_results[level] = all_results[:max_combinations]

        logger.info(f"Level {level}: Generated {len(multi_results[level])} mutations")

    return multi_results


def categorize_multi_mutations(
    multi_mutations: Dict[int, List[Dict[str, Any]]],
    regions: Optional[List[str]] = None
) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    """
    Categorize multi-mutations by region(s), mutation type, and composition.
    
    Categories:
    - single_region: All positions in the same region
    - multi_region: Positions span different regions
    - all_gatekeeper: Only gatekeeper mutations
    - all_core: Only β-core mutations
    - mixed: Combination of different mutation types
    """
    categorized = {}

    for level, mutations in multi_mutations.items():
        categorized[level] = {
            'single_region': [],
            'multi_region': [],
            'all_gatekeeper': [],
            'all_core': [],
            'mixed': [],
        }

        for item in mutations:
            is_all_gatekeeper = item.get('is_all_gatekeeper', False)
            is_all_core = item.get('is_all_core', False)
            regions_involved = item.get('regions', [])

            # Categorize by mutation type
            if is_all_gatekeeper:
                categorized[level]['all_gatekeeper'].append(item)
            elif is_all_core:
                categorized[level]['all_core'].append(item)
            else:
                categorized[level]['mixed'].append(item)

            # Categorize by region
            if not regions_involved or len(set(regions_involved)) == 1:
                categorized[level]['single_region'].append(item)
            else:
                item['region_types'] = tuple(regions_involved)
                categorized[level]['multi_region'].append(item)

    return categorized


# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

def _level_to_text(level: int) -> str:
    """Convert level number to text."""
    level_names = {
        2: 'double',
        3: 'triple',
        4: 'quadruple',
        5: 'quintuple',
        6: 'sextuple'
    }
    return level_names.get(level, f'{level}x')


def create_output_directory(
    base_path: str,
    multi_mutation_levels: Optional[List[int]] = None
) -> Path:
    """Create organized output directory structure."""
    output_dir = Path(base_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    if multi_mutation_levels:
        for level in sorted(multi_mutation_levels):
            level_dir = output_dir / f"{_level_to_text(level)}_mutations"
            level_dir.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Created directory: {level_dir}")

    return output_dir


def write_multi_mutations_by_category(
    output_dir: str,
    original_header: str,
    original_seq: str,
    categorized_mutations: Dict[int, Dict[str, List[Dict[str, Any]]]],
    include_original: bool = True
):
    """Write categorized multi-mutations to separate FASTA files."""
    output_base = Path(output_dir)

    for level in sorted(categorized_mutations.keys()):
        level_dir = output_base / f"{_level_to_text(level)}_mutations"
        level_dir.mkdir(parents=True, exist_ok=True)

        categories = categorized_mutations[level]

        # Write each category
        category_files = [
            ('single_region', 'single_region.fasta'),
            ('multi_region', 'multi_region.fasta'),
            ('all_gatekeeper', 'all_gatekeeper.fasta'),
            ('all_core', 'all_core.fasta'),
            ('mixed', 'mixed.fasta'),
        ]

        for category_name, filename in category_files:
            if categories.get(category_name):
                file_path = level_dir / filename
                mutations_list = [
                    (item['description'], item['sequence'], item.get('agg_score', 0))
                    for item in categories[category_name]
                ]
                _write_fasta_file(
                    str(file_path), original_header, original_seq,
                    mutations_list, include_original
                )
                logger.info(
                    f"Wrote {len(mutations_list)} {category_name} "
                    f"{_level_to_text(level)} mutations to {filename}"
                )


def _write_fasta_file(
    output_file: str,
    original_header: str,
    original_seq: str,
    mutations: List[Tuple[str, str, int]],
    include_original: bool = True
):
    """Write FASTA file from mutations with scores."""
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for i in range(0, len(original_seq), FASTA_LINE_LENGTH):
                f.write(f"{original_seq[i:i+FASTA_LINE_LENGTH]}\n")

        for item in mutations:
            if len(item) == 3:
                description, mutated_seq, _ = item
            else:
                description, mutated_seq = item[:2]

            protein_name = original_header[1:].strip()
            f.write(f">{protein_name}_{description}\n")
            for j in range(0, len(mutated_seq), FASTA_LINE_LENGTH):
                f.write(f"{mutated_seq[j:j+FASTA_LINE_LENGTH]}\n")


def write_fasta(
    output_file: str,
    original_header: str,
    original_seq: str,
    mutations: List[Tuple[str, str]],
    include_original: bool = True
):
    """Write results to FASTA file."""
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for i in range(0, len(original_seq), FASTA_LINE_LENGTH):
                f.write(f"{original_seq[i:i+FASTA_LINE_LENGTH]}\n")

        for description, mutated_seq in mutations:
            protein_name = original_header[1:].strip()
            f.write(f">{protein_name}_{description}\n")
            for j in range(0, len(mutated_seq), FASTA_LINE_LENGTH):
                f.write(f"{mutated_seq[j:j+FASTA_LINE_LENGTH]}\n")


# =============================================================================
# INPUT VALIDATION AND PARSING
# =============================================================================

def read_fasta(filepath: str) -> Tuple[str, str]:
    """Read a single sequence FASTA file."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        header = ""
        sequence = ""

        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    break
                header = line
            elif line and not header:
                raise ValueError("FASTA file must start with header line (>)")
            elif header:
                sequence += line.upper()

        if not header or not sequence:
            raise ValueError("Invalid FASTA file or empty sequence")

        return header, sequence

    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")
    except Exception as e:
        raise ValueError(f"Error reading FASTA file: {e}")


def validate_amino_acids(
    input_data: Union[str, List[str]],
    name: str = "amino acids",
    strict: bool = False
) -> bool:
    """Validate amino acid codes."""
    if isinstance(input_data, str):
        invalid_chars = [char for char in set(input_data.upper()) if char not in VALID_AAS]
        if invalid_chars:
            if strict:
                raise ValueError(
                    f"Invalid {name}: contains invalid characters {invalid_chars}. "
                    f"Valid amino acids: {', '.join(sorted(VALID_AAS))}"
                )
            return False
        return True

    elif isinstance(input_data, list):
        invalid = [aa for aa in input_data if len(aa) != 1 or aa.upper() not in VALID_AAS]
        if invalid:
            if strict:
                raise ValueError(
                    f"Invalid {name}: {invalid}. "
                    f"Valid amino acids: {', '.join(sorted(VALID_AAS))}"
                )
            return False
        return True

    else:
        raise TypeError(f"input_data must be str or list, got {type(input_data).__name__}")


def validate_sequence(sequence: str) -> bool:
    """Validate that the sequence contains only valid amino acid codes."""
    return validate_amino_acids(sequence, "sequence", strict=False)


def parse_region(region_str: str, seq_length: int) -> Tuple[int, int]:
    """Parse region string in format start:stop (1-indexed)."""
    try:
        if ':' not in region_str:
            raise ValueError("Region must be in format start:stop")

        start_str, stop_str = region_str.split(':')
        start = int(start_str.strip())
        stop = int(stop_str.strip())

        if start < 1 or stop > seq_length:
            raise ValueError(f"Region {start}:{stop} out of bounds (1-{seq_length})")
        if start > stop:
            raise ValueError(f"Start position {start} cannot be greater than stop position {stop}")

        return start, stop
    except ValueError as e:
        raise ValueError(f"Invalid region format '{region_str}': {e}")


def normalize_regions(
    regions: Optional[List[str]],
    seq_length: int
) -> Optional[List[str]]:
    """Normalize --regions values. Supports "all" to expand to full sequence."""
    if not regions:
        return regions

    tokens = [r.strip() for r in regions if r is not None and str(r).strip()]
    if not tokens:
        return None

    if any(t.lower() == "all" for t in tokens):
        return [f"1:{seq_length}"]

    return tokens


# =============================================================================
# HELP AND USAGE
# =============================================================================

def print_usage_example():
    """Print detailed usage example."""
    example = """
USAGE EXAMPLES:
=========================================================

1. Rule-based mutagenesis in specific regions:
   python AGGRESSOR.py protein.fasta --regions 10:20 30:40 50:60

2. Rule-based with custom mutations:
   python AGGRESSOR.py protein.fasta --regions 5:15 -m A D E

3. Rule-based with specific rules only:
   python AGGRESSOR.py protein.fasta --regions 10:30 --rules hydrophobic_aliphatic aromatic

4. Combined approach (rules + specific positions):
   python AGGRESSOR.py protein.fasta --regions 10:20 --positions 15 25 --mutations P G

5. With insertions and rule-based:
   python AGGRESSOR.py protein.fasta --regions 5:15 --insert-positions 10 --insert-aas K

6. With gatekeeping amino acids (only for edge positions):
   python AGGRESSOR.py protein.fasta --regions 10:20 --gatekeeping Y K

7. Detailed verbose output:
   python AGGRESSOR.py protein.fasta --regions 10:20 -v

8. Analyze the entire sequence (all residues):
   python AGGRESSOR.py protein.fasta --regions all

9. Generate double and triple mutations:
   python AGGRESSOR.py protein.fasta --regions 10:30 --multi-mutations 2 3

10. Multi-mutations with parallel processing:
    python AGGRESSOR.py protein.fasta --regions 10:30 --multi-mutations 2 3 --threads 4

AVAILABLE RULES:
=========================================================
• hydrophobic_aliphatic    : Triggers if ≥3 V, I, L, A, M residues within 4 positions
• aromatic                 : Triggers if ≥2 F, Y, W residues within 3 positions
• amide                    : Triggers if ≥2 Q, N residues within 3 positions
• hydrophobic_and_aromatic : Triggers if ≥2 hydrophobic-aromatic adjacent pairs 
                             OR 1 pair + at least 1 hydrophobic within 3 positions

MUTATION TYPE CLASSIFICATION:
=========================================================
• BETA_CORE  : Within identified aggregation cluster (highest risk)
• GATEKEEPER : At cluster boundary with gatekeeper AA (P, K, R, D, E)
• BOUNDARY   : At cluster boundary but not gatekeeper AA
• FLANKING   : Adjacent to APR but outside gatekeeper zone
• DIRECT     : User-specified position (no rule context)
• INSERTION  : Amino acid insertion

AGGREGATION SCORE RANKING (highest to lowest):
1. hydrophobic_aliphatic: 3
2. hydrophobic_and_aromatic: 2
3. aromatic: 2
4. amide: 1

REQUIRED PARAMETERS:
=========================================================
• input_file    : Input FASTA file (multifasta not supported)
• Either --positions OR --regions must be specified
"""
    print(example)


def print_help_info(parser: argparse.ArgumentParser):
    """Print detailed help information."""
    print("=" * 70)
    print("AGGRESSOR: RULE-BASED IN SILICO MUTAGENESIS")
    print("WITH MULTI-POINT MUTATION SUPPORT")
    print("=" * 70)
    print("\nDESCRIPTION:")
    print("Performs rule-based mutagenesis on protein sequences.")
    print("Rules apply ONLY when amino acids are clustered together.")
    print("Multiple rules can apply simultaneously to the same motif.")
    print("Supports generation of multi-point mutations (double, triple, etc.)")
    print("\n" + "=" * 70)

    print_usage_example()
    parser.print_help()

    print("\n" + "=" * 70)
    print("MULTI-POINT MUTATION FEATURES:")
    print("=" * 70)
    print("\nWhen --multi-mutations is specified, the script generates mutations at multiple")
    print("levels and organizes them in a directory structure.")
    print("\nCombinatorics control:")
    print("  • --multi-top-per-position limits variants per position (ranked by agg_score)")
    print("  • --threads enables parallel processing for large combination spaces")
    print(f"  • Maximum mutation level is {MAX_MUTATION_LEVEL}")
    print("\nOutput structure:")
    print("  mutated_sequences/")
    print("  ├── single_mutations.fasta")
    print("  ├── double_mutations/")
    print("  │   ├── single_region.fasta")
    print("  │   ├── multi_region.fasta")
    print("  │   ├── all_gatekeeper.fasta")
    print("  │   ├── all_core.fasta")
    print("  │   └── mixed.fasta")
    print("  └── triple_mutations/")
    print("      └── ...")

    print("\n" + "=" * 70)
    print("RULE DETAILS:")
    print("=" * 70)
    for rule_name, rule in RULES.items():
        print(f"\n{rule_name}:")
        print(f"  {rule['description']}")
        if rule_name == 'hydrophobic_and_aromatic':
            print(f"  Conditions:")
            print(f"    1. At least 2 hydrophobic-aromatic adjacent pairs")
            print(f"    2. OR 1 pair + at least 1 hydrophobic within 3 positions")
        else:
            print(f"  Residues: {', '.join(sorted(rule['residues']))}")
            print(f"  Min cluster size: {rule['min_cluster_size']}")
            print(f"  Max gap: {rule['max_gap']} positions")
        print(f"  Aggregation score: {rule['aggregation_score']}")

    print("\n" + "=" * 70)
    print("BIOLOGICAL REFERENCES:")
    print("=" * 70)
    print("• Rousseau et al., J Mol Biol 2006 - Gatekeeper hypothesis")
    print("• Beerten et al., FEBS Lett 2012 - APR boundary effects")
    print("• Tartaglia et al., J Mol Biol 2008 - Aggregation propensity scale")
    print("=" * 70)


def print_mutation_summary(mutations: List[Tuple[str, str]]):
    """Print summary of generated mutations."""
    if not mutations:
        print("\nNo mutations generated. Check your criteria.")
        return

    print(f"\n{'='*70}")
    print("MUTATION SUMMARY (Sorted by aggregation score)")
    print(f"{'='*70}")

    # Extract aggregation scores and group mutations
    mutation_data = []
    for desc, seq in mutations:
        agg_score = 0
        if "(agg_score=" in desc:
            try:
                agg_str = desc.split("(agg_score=")[1].split(")")[0]
                agg_score = int(agg_str)
            except (IndexError, ValueError):
                agg_score = 0
        mutation_data.append((desc, seq, agg_score))

    # Count by type
    type_counts = {mt.name: 0 for mt in MutationType}
    for desc, _, _ in mutation_data:
        for mt in MutationType:
            if mt.name in desc:
                type_counts[mt.name] += 1
                break

    print(f"\nMutation Type Breakdown:")
    for mt_name, count in type_counts.items():
        if count > 0:
            print(f"  • {mt_name}: {count}")

    print(f"\nTOTAL: {len(mutations)}")

    # Show top mutations
    if mutation_data:
        print(f"\nTOP 5 MUTATIONS BY AGGREGATION SCORE:")
        for i, (desc, _, agg_score) in enumerate(mutation_data[:5], 1):
            if len(desc) > 80:
                desc = desc[:77] + "..."
            print(f"{i}. {desc}")


def print_aggregation_summary(region_analyses: List[Dict], sequence: str):
    """Print summary of aggregation analysis results."""
    if not region_analyses:
        print("\n" + "=" * 70)
        print("AGGREGATION ANALYSIS RESULTS")
        print("=" * 70)
        print("\nNo regions analyzed. Specify regions with --regions")
        return

    print("\n" + "=" * 70)
    print("AGGREGATION ANALYSIS RESULTS")
    print("=" * 70)

    total_hotspots = 0
    total_multi_rule = 0
    all_hotspot_positions = set()

    for analysis in region_analyses:
        region_start, region_end = analysis['region']
        clusters = analysis['merged_clusters']
        multi_rule = analysis['multi_rule_clusters']
        hotspots = analysis['aggregation_hotspots']

        print(f"\nREGION {region_start}:{region_end} ({region_end - region_start + 1} residues)")
        print("-" * 70)
        print(f"Sequence: {analysis['sequence']}")
        print(f"Total clusters found: {len(clusters)}")
        print(f"Aggregation hotspot positions: {', '.join(map(str, hotspots)) if hotspots else 'None'}")

        # Show rule-by-rule breakdown
        print(f"\nRule Breakdown:")
        for rule_name, rule_data in analysis['rules'].items():
            if rule_data['condition_met']:
                num_clusters = len(rule_data['qualifying_clusters'])
                positions = rule_data['matching_positions']
                print(f"  • {rule_name}: {num_clusters} cluster(s) at positions {positions}")

        # Highlight multi-rule clusters
        if multi_rule:
            print(f"\n⚠️  HIGH AGGREGATION RISK (Multi-Rule Clusters):")
            for i, cluster in enumerate(multi_rule, 1):
                positions = cluster['positions']
                residues = ''.join(cluster['residues'])
                rules = cluster['rules']
                score = cluster['combined_aggregation_score']

                print(f"  Multi-Rule Cluster {i}:")
                print(f"    Positions: {positions}")
                print(f"    Residues:  {residues}")
                print(f"    Converging Rules: {', '.join(rules)}")
                print(f"    Combined Aggregation Score: {score}/8")

        total_hotspots += len(hotspots)
        total_multi_rule += len(multi_rule)
        all_hotspot_positions.update(hotspots)

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    print(f"Total regions analyzed: {len(region_analyses)}")
    print(f"Total aggregation hotspot positions: {total_hotspots}")
    print(f"Total multi-rule clusters (highest risk): {total_multi_rule}")

    if total_multi_rule > 0:
        print(f"\n⚠️  ATTENTION: {total_multi_rule} high-risk region(s) detected!")

    if total_hotspots > 0:
        print(f"\nRecommended mutation targets: {', '.join(map(str, sorted(all_hotspot_positions)))}")
    else:
        print("\nNo aggregation-prone hotspots detected in specified regions.")

    print("=" * 70)


# =============================================================================
# ARGUMENT PARSING AND VALIDATION
# =============================================================================

def setup_argument_parser() -> argparse.ArgumentParser:
    """Setup and return argument parser."""
    parser = argparse.ArgumentParser(
        description='Perform rule-based in silico mutagenesis on protein sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        usage='python AGGRESSOR.py <input_file> [options]'
    )

    # Required arguments
    required = parser.add_argument_group('REQUIRED ARGUMENTS')
    required.add_argument('input_file', nargs='?', help='Input FASTA file')

    # Region and rule arguments
    region_args = parser.add_argument_group('REGION-BASED MUTAGENESIS')
    region_args.add_argument(
        '-r', '--regions', type=str, nargs='+',
        help='Regions to analyze (format: start:stop) or "all"'
    )
    region_args.add_argument(
        '--rules', type=str, nargs='+',
        choices=['hydrophobic_aliphatic', 'aromatic', 'amide', 'hydrophobic_and_aromatic'],
        help='Specific rules to apply'
    )

    # Direct mutation arguments
    mutation_args = parser.add_argument_group('DIRECT MUTATIONS')
    mutation_args.add_argument(
        '-p', '--positions', type=int, nargs='+', default=[],
        help='Specific positions to mutate (1-indexed)'
    )
    mutation_args.add_argument(
        '-m', '--mutations', type=str, nargs='+', default=DEFAULT_MUTATIONS,
        help=f'Amino acids to mutate to (default: {DEFAULT_MUTATIONS})'
    )

    # Gatekeeping amino acids
    gatekeeping_args = parser.add_argument_group('GATEKEEPING AMINO ACIDS')
    gatekeeping_args.add_argument(
        '-g', '--gatekeeping', type=str, nargs='+', default=GATEKEEPING_AAS,
        help=f'Amino acids for edge positions (default: {GATEKEEPING_AAS})'
    )

    # Insertion arguments
    insertion_args = parser.add_argument_group('INSERTIONS')
    insertion_args.add_argument(
        '--insert-positions', type=int, nargs='+',
        help='Positions for insertions (before this position)'
    )
    insertion_args.add_argument(
        '--insert-aas', type=str, nargs='+',
        help='Amino acids to insert'
    )

    # Multi-point mutation arguments
    multi_mut_args = parser.add_argument_group('MULTI-POINT MUTATIONS')
    multi_mut_args.add_argument(
        '--multi-mutations', type=int, nargs='+',
        help='Levels to generate (e.g., 2 3 for double and triple)'
    )
    multi_mut_args.add_argument(
        '--multi-top-per-position', type=int, default=3,
        help='Limit variants per position (default: 3)'
    )
    multi_mut_args.add_argument(
        '--multi-output', default='mutated_sequences',
        help='Output directory for mutations (default: mutated_sequences)'
    )
    multi_mut_args.add_argument(
        '--threads', '-t', type=int, default=None,
        help='Number of parallel workers (default: CPU count - 1)'
    )

    # Aggregation analysis arguments
    agg_args = parser.add_argument_group('AGGREGATION ANALYSIS')
    agg_args.add_argument(
        '--agg-only', action='store_true',
        help='Only identify aggregation hotspots without generating mutations'
    )
    agg_args.add_argument(
        '--min-agg-score', type=int, default=4,
        help='Minimum aggregation score for hotspot (default: 4)'
    )

    # Output arguments
    output_args = parser.add_argument_group('OUTPUT')
    output_args.add_argument(
        '-o', '--output', default='mutated_sequences.fasta',
        help='Output FASTA file for single mutations'
    )
    output_args.add_argument(
        '--no-original', action='store_true',
        help='Do not include original sequence in output'
    )

    # Other options
    other_args = parser.add_argument_group('OTHER OPTIONS')
    other_args.add_argument(
        '-v', '--verbose', action='store_true',
        help='Show detailed analysis'
    )
    other_args.add_argument(
        '-h', '--help', action='store_true',
        help='Show this help message'
    )

    return parser


def validate_arguments(args) -> None:
    """Validate command line arguments."""
    if not args.agg_only and not args.positions and not args.regions and not args.insert_positions:
        print("\nERROR: You must specify at least one of:")
        print("  • --positions for direct mutations")
        print("  • --regions for rule-based mutations")
        print("  • --insert-positions for insertions")
        print("  • --agg-only for aggregation analysis only")
        print("\nUse --help for more information.")
        sys.exit(1)

    # Validate mutation list
    try:
        validate_amino_acids(args.mutations, "mutation amino acids", strict=True)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate gatekeeping amino acids
    try:
        validate_amino_acids(args.gatekeeping, "gatekeeping amino acids", strict=True)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate insertion arguments
    if bool(args.insert_positions) != bool(args.insert_aas):
        print("\nERROR: Both --insert-positions and --insert-aas must be provided together")
        sys.exit(1)

    # Validate multi-mutations
    if args.multi_mutations:
        if not isinstance(args.multi_mutations, list):
            args.multi_mutations = [args.multi_mutations]

        for level in args.multi_mutations:
            if level < 2:
                print(f"\nERROR: Multi-mutation levels must be >= 2 (got {level})")
                sys.exit(1)
            if level > MAX_MUTATION_LEVEL:
                print(f"\nERROR: Maximum mutation level is {MAX_MUTATION_LEVEL} (got {level})")
                sys.exit(1)

        if args.multi_top_per_position is not None and args.multi_top_per_position < 1:
            print("\nERROR: --multi-top-per-position must be >= 1")
            sys.exit(1)

        logger.debug(f"Multi-mutation levels: {sorted(args.multi_mutations)}")


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main entry point with comprehensive error handling."""
    parser = setup_argument_parser()
    args = parser.parse_args()

    # Setup logging BEFORE any operations
    setup_logging(verbose=args.verbose)

    if args.help or not args.input_file:
        print_help_info(parser)
        if not args.input_file:
            logger.error("Input FASTA file is required")
        sys.exit(0 if args.help else 1)

    try:
        validate_arguments(args)
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"Argument validation failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    try:
        logger.info("=" * 70)
        logger.info("AGGRESSOR: RULE-BASED MUTAGENESIS WITH OPTIMIZED CLUSTERING")
        logger.info("=" * 70)

        # Read input
        try:
            header, sequence = read_fasta(args.input_file)
        except FileNotFoundError:
            logger.error(f"Input file not found: {args.input_file}")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Invalid FASTA format: {e}")
            sys.exit(1)

        logger.info(f"Input: {args.input_file}")
        logger.info(f"Sequence length: {len(sequence)} residues")

        # Normalize regions
        try:
            args.regions = normalize_regions(args.regions, len(sequence))
        except Exception as e:
            logger.error(f"Region normalization failed: {e}")
            sys.exit(1)

        if args.regions:
            logger.info(f"Regions to analyze: {args.regions}")
        if args.positions:
            logger.info(f"Direct mutation positions: {args.positions}")
        logger.info(f"Mutations: {args.mutations}")
        logger.info(f"Gatekeeping amino acids: {args.gatekeeping}")

        if not args.agg_only:
            # Generate mutations
            mutations, region_analyses = mutate_sequence(
                sequence,
                args.positions,
                [m.upper() for m in args.mutations],
                args.regions,
                args.rules,
                args.insert_positions,
                [aa.upper() for aa in args.insert_aas] if args.insert_aas else None,
                [aa.upper() for aa in args.gatekeeping],
                args.verbose
            )

            logger.info(f"Generated {len(mutations)} mutations")

            # Determine output path
            if args.multi_mutations:
                output_base = Path(args.multi_output)
                output_base.mkdir(parents=True, exist_ok=True)

                # Write single mutations
                single_output = output_base / 'single_mutations.fasta'
                write_fasta(str(single_output), header, sequence, mutations, not args.no_original)
                logger.info(f"Single mutations written to {single_output}")

                # Generate multi-point mutations
                logger.info(f"Generating multi-point mutations (levels: {args.multi_mutations})")
                try:
                    multi_mutations = generate_multi_point_mutations(
                        mutations,
                        sequence,
                        args.multi_mutations,
                        regions=args.regions,
                        top_variants_per_position=args.multi_top_per_position,
                        n_workers=args.threads
                    )

                    # Categorize
                    categorized = categorize_multi_mutations(multi_mutations, regions=args.regions)

                    # Create output structure
                    create_output_directory(str(output_base), args.multi_mutations)

                    # Write categorized mutations
                    write_multi_mutations_by_category(
                        str(output_base),
                        header,
                        sequence,
                        categorized,
                        include_original=not args.no_original
                    )

                    logger.info(f"Multi-mutation results written to {output_base}")

                except MemoryError:
                    logger.error(
                        "Out of memory during multi-mutation generation. "
                        "Try reducing --multi-mutations levels or --multi-top-per-position"
                    )
                    sys.exit(1)

            else:
                write_fasta(args.output, header, sequence, mutations, not args.no_original)
                logger.info(f"Results written to {args.output}")

            print_mutation_summary(mutations)

        else:
            # Aggregation analysis only
            region_analyses = []
            if args.regions:
                for region_str in args.regions:
                    try:
                        start, stop = parse_region(region_str, len(sequence))
                        analysis = analyze_region(sequence, start, stop, selected_rules=args.rules)
                        region_analyses.append(analysis)
                        logger.info(f"Region {start}:{stop}: {len(analysis['merged_clusters'])} clusters found")
                    except ValueError as e:
                        logger.error(f"Error analyzing region {region_str}: {e}")

            logger.info("Aggregation analysis completed")
            print_aggregation_summary(region_analyses, sequence)

        logger.info("=" * 70)
        logger.info("ANALYSIS COMPLETE")
        logger.info("=" * 70)

    except KeyboardInterrupt:
        logger.warning("Operation cancelled by user")
        sys.exit(130)
    except MemoryError:
        logger.error("Out of memory - try reducing --multi-mutations levels or --multi-top-per-position")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

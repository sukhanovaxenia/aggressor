"""
Union-Find data structure and cluster merging operations.

Provides efficient algorithms for:
- Disjoint set union operations (Union-Find with path compression)
- Merging overlapping or proximal aggregation-prone clusters
- Cluster deduplication and grouping

The Union-Find structure enables O(α(n)) merging of overlapping clusters,
where α is the inverse Ackermann function (effectively constant time).

Architecture:
    UnionFind: Disjoint Set Union with path compression and union by rank
    merge_overlapping_clusters_unionfind: Main merging pipeline
    Cluster factory and deduplication utilities
"""
from typing import List, Dict, Set, Tuple, Union

from config import MAX_GAP_FOR_MERGING, logger
from models import Cluster, MultiRuleCluster


# =============================================================================
# UNION-FIND DATA STRUCTURE
# =============================================================================

class UnionFind:
    """
    Disjoint Set Union (Union-Find) with path compression and union by rank.

    Provides near-constant time operations for managing disjoint sets.
    Essential for efficiently merging overlapping aggregation clusters.

    Time Complexity:
        - find: O(α(n)) amortized (α = inverse Ackermann, ≤4 for any practical n)
        - union: O(α(n)) amortized
        - connected: O(α(n)) amortized

    Space Complexity: O(n) for n elements

    Usage in AGGRESSOR:
        Each identified cluster is treated as an element. When clusters
        overlap or are within gap_tolerance of each other, they are
        unioned into the same set. The final sets represent merged
        aggregation-prone regions that may span multiple overlapping
        clusters.

    Attributes:
        parent: Array tracking each element's parent (root if parent[i] == i)
        rank: Array tracking approximate tree height for union by rank
        _size: Array tracking size of each set (only valid for roots)
    """

    __slots__ = ('parent', 'rank', '_size')

    def __init__(self, n: int):
        """
        Initialize Union-Find structure for n elements.

        Each element starts in its own singleton set.

        Args:
            n: Number of elements (indexed 0 to n-1)

        Example:
            >>> uf = UnionFind(5)
            >>> uf.get_groups()
            {0: [0], 1: [1], 2: [2], 3: [3], 4: [4]}
        """
        if n < 0:
            raise ValueError(f"n must be non-negative, got {n}")

        self.parent = list(range(n))
        self.rank = [0] * n
        self._size = [1] * n

    def find(self, x: int) -> int:
        """
        Find the representative (root) of element x's set.

        Uses path compression: all nodes along the path to root
        are updated to point directly to root, flattening the tree.

        Args:
            x: Element index (0 to n-1)

        Returns:
            Root element of the set containing x

        Example:
            >>> uf = UnionFind(3)
            >>> uf.union(0, 1)
            True
            >>> uf.find(0) == uf.find(1)
            True
        """
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x: int, y: int) -> bool:
        """
        Merge the sets containing elements x and y.

        Uses union by rank: the shorter tree is attached under the
        taller tree to maintain balance and prevent degeneration
        into linked lists.

        Args:
            x: First element index
            y: Second element index

        Returns:
            True if sets were merged, False if already in same set

        Example:
            >>> uf = UnionFind(3)
            >>> uf.union(0, 1)
            True
            >>> uf.union(0, 1)
            False  # Already merged
        """
        px, py = self.find(x), self.find(y)

        if px == py:
            return False  # Already in same set

        # Union by rank: attach smaller tree under larger
        if self.rank[px] < self.rank[py]:
            px, py = py, px

        self.parent[py] = px
        self._size[px] += self._size[py]

        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1

        return True

    def connected(self, x: int, y: int) -> bool:
        """
        Check if elements x and y are in the same set.

        Args:
            x: First element index
            y: Second element index

        Returns:
            True if x and y belong to the same set

        Example:
            >>> uf = UnionFind(3)
            >>> uf.union(0, 1)
            >>> uf.connected(0, 1)
            True
            >>> uf.connected(0, 2)
            False
        """
        return self.find(x) == self.find(y)

    def set_size(self, x: int) -> int:
        """
        Get the number of elements in the set containing x.

        Args:
            x: Element index

        Returns:
            Size of the set

        Example:
            >>> uf = UnionFind(3)
            >>> uf.union(0, 1)
            >>> uf.set_size(0)
            2
        """
        return self._size[self.find(x)]

    def get_groups(self) -> Dict[int, List[int]]:
        """
        Get all disjoint sets as a dictionary.

        Returns:
            Dictionary mapping root element → list of all elements in that set

        Example:
            >>> uf = UnionFind(3)
            >>> uf.union(0, 1)
            >>> uf.get_groups()
            {0: [0, 1], 2: [2]}
        """
        groups: Dict[int, List[int]] = {}
        for i in range(len(self.parent)):
            root = self.find(i)
            if root not in groups:
                groups[root] = []
            groups[root].append(i)
        return groups


# =============================================================================
# CLUSTER FACTORY AND UTILITIES
# =============================================================================

def create_cluster(
        positions: List[int],
        sequence: str,
        rule_name: str,
        aggregation_score: int
) -> Cluster:
    """
    Factory function to create an immutable Cluster from mutable inputs.

    Handles conversion from lists to tuples and extracts residues
    from the sequence. Provides a clean interface while the Cluster
    itself remains frozen (immutable).

    Args:
        positions: List of 1-indexed positions (will be sorted)
        sequence: Full protein sequence for residue extraction
        rule_name: Name of the rule that identified this cluster
        aggregation_score: Integer score for aggregation propensity

    Returns:
        Immutable Cluster object

    Example:
        >>> seq = "ALVIF"
        >>> c = create_cluster([1, 2, 3], seq, "test", 2)
        >>> c.positions
        (1, 2, 3)
        >>> c.residues
        ('A', 'L', 'V')
    """
    sorted_positions = tuple(sorted(positions))
    residues = tuple(sequence[p - 1] for p in sorted_positions)

    return Cluster(
        positions=sorted_positions,
        residues=residues,
        rule_name=rule_name,
        aggregation_score=aggregation_score
    )


def find_position_clusters(
        positions: List[int],
        max_gap: int = 3
) -> List[List[int]]:
    """
    Group positions into clusters based on proximity.

    A simple greedy algorithm that groups consecutive positions
    within max_gap of each other. This is used for initial
    clustering within a single rule before cross-rule merging.

    Args:
        positions: Sorted list of 1-indexed positions
        max_gap: Maximum allowed gap between consecutive positions
                in the same cluster (default: 3)

    Returns:
        List of clusters, each cluster is a list of positions
        (only clusters with ≥2 positions are returned)

    Example:
        >>> find_position_clusters([1, 2, 3, 10, 11], max_gap=3)
        [[1, 2, 3], [10, 11]]
        >>> find_position_clusters([1, 5, 10], max_gap=3)
        []  # No clusters with ≥2 positions
    """
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

    # Don't forget the last cluster
    if len(current_cluster) >= 2:
        clusters.append(current_cluster)

    return clusters


# =============================================================================
# CLUSTER MERGING OPERATIONS
# =============================================================================

def merge_overlapping_clusters_unionfind(
        clusters: List[Cluster],
        sequence: str,
        gap_tolerance: int = MAX_GAP_FOR_MERGING
) -> List[Union[Cluster, MultiRuleCluster]]:
    """
    Merge overlapping or proximal clusters using Union-Find.

    Algorithm (three-phase):
    Phase 1: Build Union-Find structure from cluster overlap graph
    Phase 2: Group clusters by their Union-Find root
    Phase 3: Merge each group into a single cluster or MultiRuleCluster

    Two clusters are considered overlapping if their position ranges
    are within gap_tolerance of each other. Merged clusters that
    originated from different rules become MultiRuleCluster objects
    with boosted aggregation scores.

    Args:
        clusters: List of Cluster objects from all rules
        sequence: Full protein sequence for residue lookup
        gap_tolerance: Maximum gap between clusters to consider
                      them as overlapping (default: MAX_GAP_FOR_MERGING)

    Returns:
        List of merged clusters (Cluster or MultiRuleCluster objects).
        Non-overlapping clusters are returned unchanged.

    Example:
        >>> c1 = Cluster((1,2), ('A','L'), 'rule1', 2)
        >>> c2 = Cluster((3,4), ('V','F'), 'rule2', 1)
        >>> merged = merge_overlapping_clusters_unionfind([c1, c2], "ALVF", gap_tolerance=1)
        >>> len(merged)
        1
        >>> isinstance(merged[0], MultiRuleCluster)
        True
        >>> merged[0].rules
        ('rule1', 'rule2')
    """
    if not clusters:
        return []

    n = len(clusters)
    uf = UnionFind(n)

    # Phase 1: Identify all overlapping pairs and union them
    for i in range(n):
        for j in range(i + 1, n):
            if clusters[i].overlaps_with(clusters[j], gap_tolerance):
                uf.union(i, j)
                logger.debug(
                    f"Union clusters {i} and {j}: "
                    f"{clusters[i].rule_name} & {clusters[j].rule_name}"
                )

    # Phase 2: Group clusters by their Union-Find root
    groups = uf.get_groups()

    # Phase 3: Merge each group into a single cluster
    merged_clusters = []

    for root, indices in groups.items():
        group_clusters = [clusters[i] for i in indices]

        if len(group_clusters) == 1:
            # Single cluster, no merging needed
            merged_clusters.append(group_clusters[0])
        else:
            # Multiple overlapping clusters, merge them
            merged = _merge_cluster_group_to_object(group_clusters, sequence)
            merged_clusters.append(merged)

            logger.debug(
                f"Merged {len(group_clusters)} clusters: "
                f"{[c.rule_name for c in group_clusters]}"
            )

    return merged_clusters


def _merge_cluster_group_to_object(
        clusters: List[Cluster],
        sequence: str
) -> Union[Cluster, MultiRuleCluster]:
    """
    Merge a group of overlapping clusters into a single object.

    Combines all positions from constituent clusters, removes duplicates,
    and calculates appropriate rules and scores.

    If clusters come from different rules → MultiRuleCluster with boosted score.
    If clusters come from same rule → single Cluster with max score.

    The combined aggregation score is calculated as:
    - Multi-rule: sum of individual scores + (number of rules - 1) bonus
    - Single rule: maximum individual score

    Args:
        clusters: List of overlapping Cluster objects to merge
        sequence: Full protein sequence for residue lookup

    Returns:
        Cluster (if single rule) or MultiRuleCluster (if multiple rules)
    """
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
        # Multiple rules converged: create MultiRuleCluster with boosted score
        combined_score = total_score + len(all_rules) - 1
        return MultiRuleCluster(
            positions=sorted_positions,
            residues=residues,
            rules=tuple(sorted(all_rules)),
            combined_aggregation_score=combined_score
        )
    else:
        # Single rule: return standard Cluster
        rule_name = next(iter(all_rules))
        return Cluster(
            positions=sorted_positions,
            residues=residues,
            rule_name=rule_name,
            aggregation_score=max_score
        )
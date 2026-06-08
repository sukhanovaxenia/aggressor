"""Rule registry, default registry, and backward-compatible RULES dict."""
from typing import List, Dict, Optional

from aggressor.core.config import logger, DEFAULT_MUTATIONS
from aggressor.core.models import ClusterEvaluator, Cluster
from aggressor.rules.evaluators import (
    HydrophobicAliphaticEvaluator,
    AromaticEvaluator,
    AmideEvaluator,
    HydrophobicAromaticEvaluator,
)


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

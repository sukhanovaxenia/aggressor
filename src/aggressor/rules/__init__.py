"""Aggregation rules: evaluators, registry, and the legacy RULES dict."""
from aggressor.rules.evaluators import (
    BaseClusterEvaluator,
    HydrophobicAliphaticEvaluator,
    AromaticEvaluator,
    AmideEvaluator,
    HydrophobicAromaticEvaluator,
)
from aggressor.rules.registry import (
    RuleRegistry,
    create_default_registry,
    DEFAULT_REGISTRY,
    RULES,
)

__all__ = [
    "BaseClusterEvaluator",
    "HydrophobicAliphaticEvaluator",
    "AromaticEvaluator",
    "AmideEvaluator",
    "HydrophobicAromaticEvaluator",
    "RuleRegistry",
    "create_default_registry",
    "DEFAULT_REGISTRY",
    "RULES",
]

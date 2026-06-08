"""
AGGRESSOR: Aggregation-Guided Generation of REgion-Specific Substitution
ORiented mutations.

Rule-based in silico mutagenesis of protein sequences, targeting
aggregation-prone regions with biologically grounded gatekeeper design.

References:
    Rousseau et al., J Mol Biol 2006 (gatekeeper hypothesis)
    Beerten et al., FEBS Lett 2012 (APR boundary effects)
    Tartaglia et al., J Mol Biol 2008 (aggregation propensity scale)
"""
__version__ = "0.2.0"

from aggressor.core.models import (
    MutationType,
    Cluster,
    MultiRuleCluster,
    HydrophobicAromaticCluster,
)
from aggressor.analysis.regions import analyze_region
from aggressor.mutagenesis.engine import mutate_sequence

__all__ = [
    "__version__",
    "MutationType",
    "Cluster",
    "MultiRuleCluster",
    "HydrophobicAromaticCluster",
    "analyze_region",
    "mutate_sequence",
]

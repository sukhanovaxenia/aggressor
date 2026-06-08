"""
Configuration constants, amino acid properties, and logging setup for AGGRESSOR.

Contains all global constants used across the AGGRESSOR pipeline including
amino acid properties, aggregation propensity values, mutation parameters,
and logging configuration.

References:
- Rousseau et al., J Mol Biol 2006 (gatekeeper hypothesis)
- Beerten et al., FEBS Lett 2012 (APR boundary effects)
- Tartaglia et al., J Mol Biol 2008 (aggregation propensity scale)
"""
import logging
from typing import Dict, FrozenSet, Optional

# =============================================================================
# AMINO ACID CONSTANTS
# =============================================================================

VALID_AAS = set('ACDEFGHIKLMNPQRSTVWY')
"""Set of valid single-letter amino acid codes."""

# =============================================================================
# OUTPUT FORMATTING
# =============================================================================

FASTA_LINE_LENGTH = 60
"""Maximum line length for FASTA output format."""

# =============================================================================
# CLUSTERING PARAMETERS
# =============================================================================

MAX_GAP_FOR_MERGING = 2
"""Maximum gap between clusters to consider them for merging."""

MAX_MUTATION_LEVEL = 10
"""Hard upper limit for multi-point mutation combinations."""

# =============================================================================
# MUTATION DEFAULTS
# =============================================================================

DEFAULT_MUTATIONS = ['P', 'G', 'D', 'K']
"""Default amino acids to use for mutations (proline, glycine, aspartate, lysine)."""

# =============================================================================
# GATEKEEPING AMINO ACIDS
# =============================================================================

GATEKEEPING_AAS = ['Y']
"""
Additional gatekeeping amino acids applied only at cluster boundaries.
Tyrosine can disrupt β-sheet extension through its bulky side chain.
"""

CANONICAL_GATEKEEPER_AAS: FrozenSet[str] = frozenset({'P', 'K', 'R', 'D', 'E'})
"""
Canonical gatekeeper amino acids based on Rousseau et al., J Mol Biol 2006.
These introduce charge (K, R, D, E) or conformational constraint (P)
that disrupts β-sheet extension at aggregation-prone region boundaries.
"""

# =============================================================================
# DISTANCE THRESHOLDS FOR POSITION CLASSIFICATION
# =============================================================================

GATEKEEPER_DISTANCE = 3
"""
Maximum distance from APR boundary for gatekeeper classification.
Positions within this distance of a cluster boundary are considered
potential gatekeeper positions.
"""

FLANKING_DISTANCE = 6
"""
Maximum distance from APR boundary for flanking region classification.
Positions within this distance but beyond GATEKEEPER_DISTANCE
are classified as flanking regions.
"""

# =============================================================================
# AGGREGATION PROPENSITY SCALE
# =============================================================================

AGGREGATION_PROPENSITY: Dict[str, float] = {
    'I': 1.822,  'V': 1.594,  'L': 1.380,  'F': 1.376,
    'Y': 0.888,  'W': 0.893,  'M': 0.739,  'A': 0.411,
    'C': 0.382,  'T': 0.039,  'S': -0.228, 'G': -0.535,
    'N': -0.547, 'Q': -0.691, 'H': -0.731, 'P': -1.402,
    'D': -1.836, 'E': -1.892, 'K': -2.030, 'R': -1.814,
}
"""
Empirically-derived intrinsic aggregation propensity values for each amino acid.
Positive values indicate higher aggregation propensity.
Source: Tartaglia et al., J Mol Biol 2008.

Scale interpretation:
  > 1.0    : Strong aggregation promoters (I, V, L, F)
  0.0-1.0  : Moderate aggregation promoters (Y, W, M, A, C, T)
  -0.5-0.0 : Weak aggregation inhibitors (S, G, N)
  < -1.0   : Strong aggregation inhibitors (P, D, E, K, R)
"""

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

# Module-level logger (used within this module only)
logger = logging.getLogger(__name__)

# Default logging format strings
CONSOLE_LOG_FORMAT = '%(asctime)s | %(levelname)-8s | %(message)s'
CONSOLE_DATE_FORMAT = '%H:%M:%S'

FILE_LOG_FORMAT = '%(asctime)s | %(levelname)-8s | %(funcName)s | %(message)s'
FILE_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'


def setup_logging(
    verbose: bool = False,
    log_file: Optional[str] = None
) -> None:
    """
    Configure logging for the AGGRESSOR mutagenesis pipeline.

    Sets up console and optional file logging with appropriate formatting.
    Console handler shows INFO by default (DEBUG if verbose=True).
    File handler always logs at DEBUG level for detailed traceability.

    Args:
        verbose: If True, set console logging to DEBUG level
        log_file: Optional path to a file for detailed debug logs

    Example:
        >>> setup_logging(verbose=True, log_file="aggressor.log")
        >>> logger = logging.getLogger(__name__)
        >>> logger.info("Pipeline started")
    """
    log_level = logging.DEBUG if verbose else logging.INFO

    # Console formatter - compact for readability
    console_formatter = logging.Formatter(
        CONSOLE_LOG_FORMAT,
        datefmt=CONSOLE_DATE_FORMAT
    )

    # File formatter - includes function name for debugging
    file_formatter = logging.Formatter(
        FILE_LOG_FORMAT,
        datefmt=FILE_DATE_FORMAT
    )

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Clear existing handlers to avoid duplicates on repeated calls
    root_logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)

    # File handler (optional, always DEBUG level)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)
        logger.debug(f"Logging to file: {log_file}")

    logger.debug(f"Logging configured (verbose={verbose}, file={log_file is not None})")
"""
Command-line argument parsing, validation, and help display for AGGRESSOR.
"""
import argparse
import sys
from typing import List, Union

from aggressor.core.config import (
    GATEKEEPING_AAS,
    DEFAULT_MUTATIONS,
    MAX_MUTATION_LEVEL,
    VALID_AAS,
)
from aggressor.rules import RULES


# =============================================================================
# ARGUMENT PARSER
# =============================================================================

def setup_argument_parser() -> argparse.ArgumentParser:
    """Configure and return the AGGRESSOR argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            'AGGRESSOR: Aggregation-Guided Generation of REgion-Specific '
            'Substitution ORiented mutations. Rule-based in silico '
            'mutagenesis on protein sequences with multi-point support.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        usage='aggressor <input_file> [options]',
    )

    required = parser.add_argument_group('REQUIRED ARGUMENTS')
    required.add_argument('input_file', nargs='?', help='Input FASTA file')

    region_args = parser.add_argument_group('REGION-BASED MUTAGENESIS')
    region_args.add_argument(
        '-r', '--regions', type=str, nargs='+',
        help='Regions to analyze (format: start:stop) or "all"',
    )
    region_args.add_argument(
        '--rules', type=str, nargs='+',
        choices=['hydrophobic_aliphatic', 'aromatic', 'amide', 'hydrophobic_and_aromatic'],
        help='Specific rules to apply',
    )

    mutation_args = parser.add_argument_group('DIRECT MUTATIONS')
    mutation_args.add_argument(
        '-p', '--positions', type=int, nargs='+', default=[],
        help='Specific positions to mutate (1-indexed)',
    )
    mutation_args.add_argument(
        '-m', '--mutations', type=str, nargs='+', default=DEFAULT_MUTATIONS,
        help=f'Amino acids to mutate to (default: {DEFAULT_MUTATIONS})',
    )

    gatekeeping_args = parser.add_argument_group('GATEKEEPING AMINO ACIDS')
    gatekeeping_args.add_argument(
        '-g', '--gatekeeping', type=str, nargs='+', default=GATEKEEPING_AAS,
        help=f'Amino acids for APR boundary positions (default: {GATEKEEPING_AAS})',
    )
    gatekeeping_args.add_argument(
        '--max-gatekeepers-per-apr', type=int, default=None,
        help=(
            'Optional cap on the number of gatekept boundary slots per APR. '
            'By default the count is chosen automatically from APR geometry '
            '(both flanks). Use this only to limit variants for a synthesis '
            'budget; slots are then ranked by predicted propensity reduction.'
        ),
    )

    insertion_args = parser.add_argument_group('INSERTIONS')
    insertion_args.add_argument(
        '--insert-positions', type=int, nargs='+',
        help='Positions for insertions (before this position)',
    )
    insertion_args.add_argument(
        '--insert-aas', type=str, nargs='+', help='Amino acids to insert',
    )

    multi_args = parser.add_argument_group('MULTI-POINT MUTATIONS')
    multi_args.add_argument(
        '--multi-mutations', type=int, nargs='+',
        help='Levels to generate (e.g., 2 3 for double and triple)',
    )
    multi_args.add_argument(
        '--multi-top-per-position', type=int, default=3,
        help='Limit variants per position (default: 3)',
    )
    multi_args.add_argument(
        '--multi-output', default='mutated_sequences',
        help='Output directory for multi-mutations (default: mutated_sequences)',
    )
    multi_args.add_argument(
        '--threads', '-t', type=int, default=None,
        help='Number of parallel workers (default: CPU count - 1)',
    )

    agg_args = parser.add_argument_group('AGGREGATION ANALYSIS')
    agg_args.add_argument(
        '--agg-only', action='store_true',
        help='Only identify aggregation hotspots without generating mutations',
    )

    output_args = parser.add_argument_group('OUTPUT')
    output_args.add_argument(
        '-o', '--output', default='mutated_sequences.fasta',
        help='Output FASTA file for single mutations',
    )
    output_args.add_argument(
        '--no-original', action='store_true',
        help='Do not include the original sequence in output',
    )

    other_args = parser.add_argument_group('OTHER OPTIONS')
    other_args.add_argument(
        '-v', '--verbose', action='store_true', help='Show detailed analysis',
    )
    other_args.add_argument(
        '--version', action='store_true', help='Show version and exit',
    )
    other_args.add_argument(
        '-h', '--help', action='store_true', help='Show this help message',
    )
    return parser


# =============================================================================
# VALIDATION
# =============================================================================

def validate_amino_acids(
        input_data: Union[str, List[str]],
        name: str = "amino acids",
        strict: bool = False,
) -> bool:
    """Validate amino-acid codes against the canonical 20-letter alphabet."""
    if isinstance(input_data, str):
        invalid = [c for c in set(input_data.upper()) if c not in VALID_AAS]
    elif isinstance(input_data, list):
        invalid = [aa for aa in input_data if len(aa) != 1 or aa.upper() not in VALID_AAS]
    else:
        raise TypeError(f"input_data must be str or list, got {type(input_data).__name__}")

    if invalid:
        if strict:
            raise ValueError(
                f"Invalid {name}: {invalid}. "
                f"Valid amino acids: {', '.join(sorted(VALID_AAS))}"
            )
        return False
    return True


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate parsed arguments; exit non-zero on user error."""
    if (not args.agg_only and not args.positions
            and not args.regions and not args.insert_positions):
        print("\nERROR: You must specify at least one of:", file=sys.stderr)
        print("  - --positions for direct mutations", file=sys.stderr)
        print("  - --regions for rule-based mutations", file=sys.stderr)
        print("  - --insert-positions for insertions", file=sys.stderr)
        print("  - --agg-only for aggregation analysis only", file=sys.stderr)
        sys.exit(1)

    try:
        validate_amino_acids(args.mutations, "mutation amino acids", strict=True)
        validate_amino_acids(args.gatekeeping, "gatekeeping amino acids", strict=True)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if bool(args.insert_positions) != bool(args.insert_aas):
        print("\nERROR: --insert-positions and --insert-aas must be given together",
              file=sys.stderr)
        sys.exit(1)

    if args.max_gatekeepers_per_apr is not None and args.max_gatekeepers_per_apr < 0:
        print("\nERROR: --max-gatekeepers-per-apr must be >= 0", file=sys.stderr)
        sys.exit(1)

    if args.multi_mutations:
        for level in args.multi_mutations:
            if level < 2:
                print(f"\nERROR: Multi-mutation levels must be >= 2 (got {level})",
                      file=sys.stderr)
                sys.exit(1)
            if level > MAX_MUTATION_LEVEL:
                print(f"\nERROR: Maximum mutation level is {MAX_MUTATION_LEVEL} (got {level})",
                      file=sys.stderr)
                sys.exit(1)
        if args.multi_top_per_position is not None and args.multi_top_per_position < 1:
            print("\nERROR: --multi-top-per-position must be >= 1", file=sys.stderr)
            sys.exit(1)


# =============================================================================
# HELP DISPLAY
# =============================================================================

def print_help_info(parser: argparse.ArgumentParser) -> None:
    """Print detailed help including examples, rules, and references."""
    print("=" * 70)
    print("AGGRESSOR: AGGREGATION-GUIDED IN SILICO MUTAGENESIS")
    print("=" * 70)
    print("\nRules apply ONLY when aggregation-prone residues are clustered.")
    print("Multiple rules may converge on the same motif (multi-rule clusters).")
    print("Supports single and multi-point (double, triple, ...) mutations.\n")
    parser.print_help()

    print("\n" + "=" * 70)
    print("RULE DETAILS")
    print("=" * 70)
    for rule_name, rule in RULES.items():
        print(f"\n{rule_name}:")
        print(f"  {rule['description']}")
        if rule_name == 'hydrophobic_and_aromatic':
            print("  Condition A: >=2 hydrophobic-aromatic adjacent pairs")
            print("  Condition B: 1 pair + >=1 hydrophobic within 3 positions")
        else:
            print(f"  Residues: {', '.join(sorted(rule['residues']))}")
            print(f"  Min cluster size: {rule['min_cluster_size']}")
            print(f"  Max gap: {rule['max_gap']} positions")
        print(f"  Aggregation score: {rule['aggregation_score']}")

    print("\n" + "=" * 70)
    print("MUTATION TYPE CLASSIFICATION")
    print("=" * 70)
    print("  BETA_CORE  : within an identified aggregation cluster")
    print("  GATEKEEPER : at/near an APR boundary, mutated to P/K/R/D/E")
    print("  BOUNDARY   : at an APR boundary but not to a gatekeeper AA")
    print("  FLANKING   : adjacent to an APR, outside the gatekeeper zone")
    print("  DIRECT     : user-specified position (no rule context)")
    print("  INSERTION  : amino-acid insertion")

    print("\n" + "=" * 70)
    print("REFERENCES")
    print("=" * 70)
    print("  Rousseau et al., J Mol Biol 2006   - gatekeeper hypothesis")
    print("  Beerten et al., FEBS Lett 2012     - APR boundary effects")
    print("  Tartaglia et al., J Mol Biol 2008  - aggregation propensity scale")
    print("=" * 70)

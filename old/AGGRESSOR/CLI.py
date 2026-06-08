"""
Command-line interface for AGGRESSOR mutagenesis pipeline.

Provides argument parsing, input validation, help display,
and output formatting. Acts as the main entry point orchestrating
all pipeline components.

Architecture:
    setup_argument_parser() — configure argparse
    validate_arguments() — input validation
    print_help_info() — detailed help display
    print_mutation_summary() — results summary
    print_aggregation_summary() — analysis summary
    main() — entry point

Usage:
    python AGGRESSOR.py <input_file> [options]
    python -m aggressor.CLI <input_file> [options]
"""
import argparse
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Optional

from config import (
    GATEKEEPING_AAS,
    FASTA_LINE_LENGTH,
    DEFAULT_MUTATIONS,
    MAX_MUTATION_LEVEL,
    setup_logging,
    logger
)
from models import MutationType
from rules import RULES
from analysis import analyze_region
from mutagenesis import mutate_sequence, generate_multi_point_mutations, categorize_multi_mutations


# =============================================================================
# FASTA I/O
# =============================================================================

def read_fasta(filepath: str) -> Tuple[str, str]:
    """
    Read a single-sequence FASTA file.

    Args:
        filepath: Path to FASTA file

    Returns:
        Tuple of (header_line, sequence)

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid or empty

    Example:
        >>> header, seq = read_fasta("protein.fasta")
        >>> header
        '>sp|P04637|P53_HUMAN'
        >>> len(seq)
        393
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        header = ""
        sequence = ""

        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    # Already have header — single sequence mode
                    break
                header = line
            elif line and not header:
                raise ValueError(
                    "FASTA file must start with header line (>)"
                )
            elif header:
                sequence += line.upper()

        if not header or not sequence:
            raise ValueError("Invalid FASTA file or empty sequence")

        return header, sequence

    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")
    except Exception as e:
        raise ValueError(f"Error reading FASTA file: {e}")


def write_fasta(
        output_file: str,
        original_header: str,
        original_seq: str,
        mutations: List[Tuple[str, str]],
        include_original: bool = True
) -> None:
    """
    Write mutations to FASTA format file.

    Args:
        output_file: Path to output FASTA file
        original_header: Original FASTA header line
        original_seq: Original protein sequence
        mutations: List of (description, sequence) tuples
        include_original: Whether to include original sequence first
    """
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for i in range(0, len(original_seq), FASTA_LINE_LENGTH):
                f.write(f"{original_seq[i:i + FASTA_LINE_LENGTH]}\n")

        protein_name = original_header[1:].strip()
        for description, mutated_seq in mutations:
            f.write(f">{protein_name}_{description}\n")
            for j in range(0, len(mutated_seq), FASTA_LINE_LENGTH):
                f.write(f"{mutated_seq[j:j + FASTA_LINE_LENGTH]}\n")


def _write_fasta_file(
        output_file: str,
        original_header: str,
        original_seq: str,
        mutations: List[Tuple[str, str, int]],
        include_original: bool = True
) -> None:
    """
    Write mutations with scores to FASTA file.

    Args:
        output_file: Path to output FASTA file
        original_header: Original FASTA header line
        original_seq: Original protein sequence
        mutations: List of (description, sequence, score) tuples
        include_original: Whether to include original sequence
    """
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for i in range(0, len(original_seq), FASTA_LINE_LENGTH):
                f.write(f"{original_seq[i:i + FASTA_LINE_LENGTH]}\n")

        protein_name = original_header[1:].strip()
        for item in mutations:
            description = item[0]
            mutated_seq = item[1]
            f.write(f">{protein_name}_{description}\n")
            for j in range(0, len(mutated_seq), FASTA_LINE_LENGTH):
                f.write(f"{mutated_seq[j:j + FASTA_LINE_LENGTH]}\n")


# =============================================================================
# ARGUMENT PARSER
# =============================================================================

def setup_argument_parser() -> argparse.ArgumentParser:
    """
    Configure and return the argument parser.

    Returns:
        Configured ArgumentParser with all option groups
    """
    parser = argparse.ArgumentParser(
        description=(
            'AGGRESSOR: Aggregation-Guided Generation of REgion-Specific '
            'Substitution ORiented mutations.\n'
            'Performs rule-based in silico mutagenesis on protein sequences '
            'with multi-point mutation support.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        usage='python AGGRESSOR.py <input_file> [options]'
    )

    # Required arguments
    required = parser.add_argument_group('REQUIRED ARGUMENTS')
    required.add_argument(
        'input_file', nargs='?',
        help='Input FASTA file containing protein sequence'
    )

    # Region and rule arguments
    region_args = parser.add_argument_group('REGION-BASED MUTAGENESIS')
    region_args.add_argument(
        '-r', '--regions', type=str, nargs='+',
        help='Regions to analyze (format: start:stop) or "all"'
    )
    region_args.add_argument(
        '--rules', type=str, nargs='+',
        choices=[
            'hydrophobic_aliphatic', 'aromatic',
            'amide', 'hydrophobic_and_aromatic'
        ],
        help='Specific rules to apply'
    )

    # Direct mutation arguments
    mutation_args = parser.add_argument_group('DIRECT MUTATIONS')
    mutation_args.add_argument(
        '-p', '--positions', type=int, nargs='+', default=[],
        help='Specific positions to mutate (1-indexed)'
    )
    mutation_args.add_argument(
        '-m', '--mutations', type=str, nargs='+',
        default=DEFAULT_MUTATIONS,
        help=f'Amino acids to mutate to (default: {DEFAULT_MUTATIONS})'
    )

    # Gatekeeping amino acids
    gatekeeping_args = parser.add_argument_group('GATEKEEPING AMINO ACIDS')
    gatekeeping_args.add_argument(
        '-g', '--gatekeeping', type=str, nargs='+',
        default=GATEKEEPING_AAS,
        help=(
            f'Amino acids for edge positions '
            f'(default: {GATEKEEPING_AAS})'
        )
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
    multi_args = parser.add_argument_group('MULTI-POINT MUTATIONS')
    multi_args.add_argument(
        '--multi-mutations', type=int, nargs='+',
        help='Levels to generate (e.g., 2 3 for double and triple)'
    )
    multi_args.add_argument(
        '--multi-top-per-position', type=int, default=3,
        help='Limit variants per position (default: 3)'
    )
    multi_args.add_argument(
        '--multi-output', default='mutated_sequences',
        help='Output directory for mutations (default: mutated_sequences)'
    )
    multi_args.add_argument(
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


# =============================================================================
# HELP AND USAGE DISPLAY
# =============================================================================

def print_usage_example() -> None:
    """Print detailed usage examples."""
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


def print_help_info(parser: argparse.ArgumentParser) -> None:
    """
    Print detailed help information including examples and rule descriptions.

    Args:
        parser: Configured ArgumentParser
    """
    print("=" * 70)
    print("AGGRESSOR: AGGREGATION-GUIDED IN SILICO MUTAGENESIS")
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


# =============================================================================
# INPUT VALIDATION
# =============================================================================

def validate_amino_acids(
        input_data: str | List[str],
        name: str = "amino acids",
        strict: bool = False
) -> bool:
    """
    Validate amino acid codes against valid set.

    Args:
        input_data: String of amino acids or list of amino acid codes
        name: Description for error messages
        strict: If True, raise ValueError on invalid input

    Returns:
        True if valid, False otherwise (when strict=False)

    Raises:
        ValueError: If strict=True and invalid amino acids found
        TypeError: If input_data is not str or list
    """
    from config import VALID_AAS

    if isinstance(input_data, str):
        invalid_chars = [
            char for char in set(input_data.upper())
            if char not in VALID_AAS
        ]
        if invalid_chars:
            if strict:
                raise ValueError(
                    f"Invalid {name}: contains invalid characters "
                    f"{invalid_chars}. "
                    f"Valid amino acids: {', '.join(sorted(VALID_AAS))}"
                )
            return False
        return True

    elif isinstance(input_data, list):
        invalid = [
            aa for aa in input_data
            if len(aa) != 1 or aa.upper() not in VALID_AAS
        ]
        if invalid:
            if strict:
                raise ValueError(
                    f"Invalid {name}: {invalid}. "
                    f"Valid amino acids: {', '.join(sorted(VALID_AAS))}"
                )
            return False
        return True

    else:
        raise TypeError(
            f"input_data must be str or list, got {type(input_data).__name__}"
        )


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate command line arguments.

    Args:
        args: Parsed command line arguments

    Raises:
        SystemExit: If validation fails
    """
    if (
            not args.agg_only
            and not args.positions
            and not args.regions
            and not args.insert_positions
    ):
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

    # Validate insertion arguments come in pairs
    if bool(args.insert_positions) != bool(args.insert_aas):
        print(
            "\nERROR: Both --insert-positions and --insert-aas "
            "must be provided together",
            file=sys.stderr
        )
        sys.exit(1)

    # Validate multi-mutation levels
    if args.multi_mutations:
        for level in args.multi_mutations:
            if level < 2:
                print(
                    f"\nERROR: Multi-mutation levels must be >= 2 (got {level})",
                    file=sys.stderr
                )
                sys.exit(1)
            if level > MAX_MUTATION_LEVEL:
                print(
                    f"\nERROR: Maximum mutation level is "
                    f"{MAX_MUTATION_LEVEL} (got {level})",
                    file=sys.stderr
                )
                sys.exit(1)

        if args.multi_top_per_position is not None and args.multi_top_per_position < 1:
            print(
                "\nERROR: --multi-top-per-position must be >= 1",
                file=sys.stderr
            )
            sys.exit(1)


def parse_region(region_str: str, seq_length: int) -> Tuple[int, int]:
    """
    Parse region string in format start:stop (1-indexed).

    Args:
        region_str: Region specification (e.g., "10:50")
        seq_length: Length of the full sequence

    Returns:
        Tuple of (start, stop) positions

    Raises:
        ValueError: If format is invalid or out of bounds
    """
    try:
        if ':' not in region_str:
            raise ValueError("Region must be in format start:stop")

        start_str, stop_str = region_str.split(':')
        start = int(start_str.strip())
        stop = int(stop_str.strip())

        if start < 1 or stop > seq_length:
            raise ValueError(
                f"Region {start}:{stop} out of bounds (1-{seq_length})"
            )
        if start > stop:
            raise ValueError(
                f"Start position {start} cannot be greater than "
                f"stop position {stop}"
            )

        return start, stop
    except ValueError as e:
        raise ValueError(f"Invalid region format '{region_str}': {e}")


def normalize_regions(
        regions: Optional[List[str]],
        seq_length: int
) -> Optional[List[str]]:
    """
    Normalize --regions values. Supports "all" to expand to full sequence.

    Args:
        regions: List of region strings or None
        seq_length: Length of the full sequence

    Returns:
        Normalized list of region strings or None
    """
    if not regions:
        return None

    tokens = [
        r.strip() for r in regions
        if r is not None and str(r).strip()
    ]
    if not tokens:
        return None

    if any(t.lower() == "all" for t in tokens):
        return [f"1:{seq_length}"]

    return tokens


# =============================================================================
# OUTPUT AND SUMMARIES
# =============================================================================

def print_mutation_summary(mutations: List[Tuple[str, str]]) -> None:
    """
    Print summary of generated mutations with statistics.

    Args:
        mutations: List of (description, sequence) tuples
    """
    if not mutations:
        print("\nNo mutations generated. Check your criteria.")
        return

    print(f"\n{'=' * 70}")
    print("MUTATION SUMMARY (Sorted by aggregation score)")
    print(f"{'=' * 70}")

    # Extract aggregation scores
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

    # Count by mutation type
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


def print_aggregation_summary(
        region_analyses: List[Dict],
        sequence: str
) -> None:
    """
    Print summary of aggregation analysis results.

    Args:
        region_analyses: List of region analysis dictionaries
        sequence: Full protein sequence
    """
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

        print(
            f"\nREGION {region_start}:{region_end} "
            f"({region_end - region_start + 1} residues)"
        )
        print("-" * 70)
        print(f"Sequence: {analysis['sequence']}")
        print(f"Total clusters found: {len(clusters)}")
        print(
            f"Aggregation hotspot positions: "
            f"{', '.join(map(str, hotspots)) if hotspots else 'None'}"
        )

        # Rule-by-rule breakdown
        print(f"\nRule Breakdown:")
        for rule_name, rule_data in analysis['rules'].items():
            if rule_data['condition_met']:
                num_clusters = len(rule_data['qualifying_clusters'])
                positions = rule_data['matching_positions']
                print(
                    f"  • {rule_name}: {num_clusters} cluster(s) "
                    f"at positions {positions}"
                )

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
        print(
            f"\nRecommended mutation targets: "
            f"{', '.join(map(str, sorted(all_hotspot_positions)))}"
        )
    else:
        print("\nNo aggregation-prone hotspots detected in specified regions.")

    print("=" * 70)


# =============================================================================
# MULTI-MUTATION OUTPUT
# =============================================================================

def _level_to_text(level: int) -> str:
    """
    Convert mutation level number to text.

    Args:
        level: Number of simultaneous mutations

    Returns:
        Text description (e.g., 2 → "double", 3 → "triple")
    """
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
    """
    Create organized output directory structure.

    Args:
        base_path: Base directory path
        multi_mutation_levels: List of mutation levels

    Returns:
        Path object for the output directory
    """
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
        categorized_mutations: Dict[int, Dict[str, List[Dict]]],
        include_original: bool = True
) -> None:
    """
    Write categorized multi-mutations to separate FASTA files.

    Args:
        output_dir: Base output directory
        original_header: Original FASTA header
        original_seq: Original sequence
        categorized_mutations: Categorized mutations by level and type
        include_original: Include original sequence in output
    """
    output_base = Path(output_dir)

    for level in sorted(categorized_mutations.keys()):
        level_dir = output_base / f"{_level_to_text(level)}_mutations"
        level_dir.mkdir(parents=True, exist_ok=True)

        categories = categorized_mutations[level]

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
                    (
                        item['description'],
                        item['sequence'],
                        item.get('agg_score', 0)
                    )
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


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main() -> None:
    """Main entry point for AGGRESSOR pipeline."""
    parser = setup_argument_parser()
    args = parser.parse_args()

    # Setup logging before any operations
    setup_logging(verbose=args.verbose)

    if args.help or not args.input_file:
        print_help_info(parser)
        if not args.input_file:
            logger.error("Input FASTA file is required")
            sys.exit(1)
        sys.exit(0)

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
            # Generate single mutations
            mutations, region_analyses = mutate_sequence(
                sequence,
                args.positions,
                [m.upper() for m in args.mutations],
                args.regions,
                args.rules,
                args.insert_positions,
                (
                    [aa.upper() for aa in args.insert_aas]
                    if args.insert_aas else None
                ),
                [aa.upper() for aa in args.gatekeeping],
                args.verbose
            )

            logger.info(f"Generated {len(mutations)} mutations")

            # Handle multi-point mutations
            if args.multi_mutations:
                output_base = Path(args.multi_output)
                output_base.mkdir(parents=True, exist_ok=True)

                # Write single mutations
                single_output = output_base / 'single_mutations.fasta'
                write_fasta(
                    str(single_output), header, sequence,
                    mutations, not args.no_original
                )
                logger.info(f"Single mutations written to {single_output}")

                # Generate multi-point mutations
                logger.info(
                    f"Generating multi-point mutations "
                    f"(levels: {args.multi_mutations})"
                )
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
                    categorized = categorize_multi_mutations(
                        multi_mutations, regions=args.regions
                    )

                    # Create output structure
                    create_output_directory(
                        str(output_base), args.multi_mutations
                    )

                    # Write categorized mutations
                    write_multi_mutations_by_category(
                        str(output_base), header, sequence,
                        categorized, include_original=not args.no_original
                    )

                    logger.info(
                        f"Multi-mutation results written to {output_base}"
                    )

                except MemoryError:
                    logger.error(
                        "Out of memory during multi-mutation generation. "
                        "Try reducing --multi-mutations levels or "
                        "--multi-top-per-position"
                    )
                    sys.exit(1)

            else:
                # Single mutations only
                write_fasta(
                    args.output, header, sequence,
                    mutations, not args.no_original
                )
                logger.info(f"Results written to {args.output}")

            print_mutation_summary(mutations)

        else:
            # Aggregation analysis only
            region_analyses = []
            if args.regions:
                for region_str in args.regions:
                    try:
                        start, stop = parse_region(
                            region_str, len(sequence)
                        )
                        analysis = analyze_region(
                            sequence, start, stop,
                            selected_rules=args.rules
                        )
                        region_analyses.append(analysis)
                        logger.info(
                            f"Region {start}:{stop}: "
                            f"{len(analysis['merged_clusters'])} clusters found"
                        )
                    except ValueError as e:
                        logger.error(
                            f"Error analyzing region {region_str}: {e}"
                        )

            logger.info("Aggregation analysis completed")
            print_aggregation_summary(region_analyses, sequence)

        logger.info("=" * 70)
        logger.info("ANALYSIS COMPLETE")
        logger.info("=" * 70)

    except KeyboardInterrupt:
        logger.warning("Operation cancelled by user")
        sys.exit(130)
    except MemoryError:
        logger.error(
            "Out of memory - try reducing --multi-mutations levels "
            "or --multi-top-per-position"
        )
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
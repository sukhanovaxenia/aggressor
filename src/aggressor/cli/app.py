"""
AGGRESSOR command-line application: orchestrates the pipeline and prints
human-readable summaries.
"""
import sys
from pathlib import Path
from typing import List, Dict

from aggressor import __version__
from aggressor.core.config import setup_logging, logger
from aggressor.core.models import MutationType
from aggressor.analysis.regions import analyze_region
from aggressor.mutagenesis.engine import (
    mutate_sequence,
    generate_multi_point_mutations,
    categorize_multi_mutations,
)
from aggressor.io.fasta import (
    read_fasta,
    write_fasta,
    parse_region,
    normalize_regions,
    create_output_directory,
    write_multi_mutations_by_category,
)
from aggressor.cli.parser import (
    setup_argument_parser,
    validate_arguments,
    print_help_info,
)


# =============================================================================
# SUMMARIES
# =============================================================================

def print_mutation_summary(mutations: List[tuple]) -> None:
    """Print a breakdown of generated single mutations by type."""
    if not mutations:
        print("\nNo mutations generated. Check your criteria.")
        return

    print("\n" + "=" * 70)
    print("MUTATION SUMMARY (sorted by aggregation score)")
    print("=" * 70)

    type_counts = {mt.name: 0 for mt in MutationType}
    for desc, _ in mutations:
        for mt in MutationType:
            if mt.name in desc:
                type_counts[mt.name] += 1
                break

    print("\nMutation type breakdown:")
    for name, count in type_counts.items():
        if count:
            print(f"  - {name}: {count}")
    print(f"\nTOTAL: {len(mutations)}")

    print("\nTop 5 mutations by aggregation score:")
    for i, (desc, _) in enumerate(mutations[:5], 1):
        shown = desc if len(desc) <= 80 else desc[:77] + "..."
        print(f"  {i}. {shown}")


def print_aggregation_summary(region_analyses: List[Dict], sequence: str) -> None:
    """Print aggregation-analysis results for --agg-only mode."""
    print("\n" + "=" * 70)
    print("AGGREGATION ANALYSIS RESULTS")
    print("=" * 70)
    if not region_analyses:
        print("\nNo regions analyzed. Specify regions with --regions")
        return

    total_hotspots = 0
    total_multi_rule = 0
    all_hotspots = set()

    for analysis in region_analyses:
        rstart, rend = analysis['region']
        clusters = analysis['merged_clusters']
        multi_rule = analysis['multi_rule_clusters']
        hotspots = analysis['aggregation_hotspots']

        print(f"\nREGION {rstart}:{rend} ({rend - rstart + 1} residues)")
        print("-" * 70)
        print(f"Sequence: {analysis['sequence']}")
        print(f"Total clusters found: {len(clusters)}")
        print(f"Hotspot positions: {', '.join(map(str, hotspots)) if hotspots else 'None'}")

        print("\nRule breakdown:")
        for rule_name, rule_data in analysis['rules'].items():
            if rule_data['condition_met']:
                n = len(rule_data['qualifying_clusters'])
                print(f"  - {rule_name}: {n} cluster(s) at {rule_data['matching_positions']}")

        if multi_rule:
            print("\n[!] HIGH AGGREGATION RISK (multi-rule clusters):")
            for i, c in enumerate(multi_rule, 1):
                print(f"  Multi-rule cluster {i}: positions {c['positions']}")
                print(f"    Residues: {''.join(c['residues'])}")
                print(f"    Converging rules: {', '.join(c['rules'])}")
                print(f"    Combined score: {c['combined_aggregation_score']}")

        total_hotspots += len(hotspots)
        total_multi_rule += len(multi_rule)
        all_hotspots.update(hotspots)

    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    print(f"Regions analyzed: {len(region_analyses)}")
    print(f"Hotspot positions: {total_hotspots}")
    print(f"Multi-rule clusters (highest risk): {total_multi_rule}")
    if all_hotspots:
        print(f"\nRecommended targets: {', '.join(map(str, sorted(all_hotspots)))}")
    print("=" * 70)


# =============================================================================
# MAIN
# =============================================================================

def main(argv: List[str] = None) -> None:
    """Entry point for the ``aggressor`` console script and ``python -m aggressor``."""
    parser = setup_argument_parser()
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        print(f"AGGRESSOR {__version__}")
        sys.exit(0)

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
    except Exception as e:  # pragma: no cover - defensive
        logger.error(f"Argument validation failed: {e}")
        sys.exit(1)

    try:
        logger.info("=" * 70)
        logger.info("AGGRESSOR: RULE-BASED MUTAGENESIS")
        logger.info("=" * 70)

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

        args.regions = normalize_regions(args.regions, len(sequence))
        if args.regions:
            logger.info(f"Regions to analyze: {args.regions}")
        if args.positions:
            logger.info(f"Direct mutation positions: {args.positions}")
        logger.info(f"Mutations: {args.mutations}")
        logger.info(f"Gatekeeping amino acids: {args.gatekeeping}")
        if args.max_gatekeepers_per_apr is not None:
            logger.info(f"Max gatekeepers per APR: {args.max_gatekeepers_per_apr}")

        if not args.agg_only:
            mutations, _ = mutate_sequence(
                sequence,
                args.positions,
                [m.upper() for m in args.mutations],
                args.regions,
                args.rules,
                args.insert_positions,
                [aa.upper() for aa in args.insert_aas] if args.insert_aas else None,
                [aa.upper() for aa in args.gatekeeping],
                args.verbose,
                max_gatekeepers_per_apr=args.max_gatekeepers_per_apr,
            )
            logger.info(f"Generated {len(mutations)} mutations")

            if args.multi_mutations:
                output_base = Path(args.multi_output)
                output_base.mkdir(parents=True, exist_ok=True)
                single_output = output_base / 'single_mutations.fasta'
                write_fasta(str(single_output), header, sequence,
                            mutations, not args.no_original)
                logger.info(f"Single mutations written to {single_output}")

                logger.info(f"Generating multi-point mutations (levels: {args.multi_mutations})")
                try:
                    multi = generate_multi_point_mutations(
                        mutations, sequence, args.multi_mutations,
                        regions=args.regions,
                        top_variants_per_position=args.multi_top_per_position,
                        n_workers=args.threads,
                    )
                    categorized = categorize_multi_mutations(multi, regions=args.regions)
                    create_output_directory(str(output_base), args.multi_mutations)
                    write_multi_mutations_by_category(
                        str(output_base), header, sequence,
                        categorized, include_original=not args.no_original,
                    )
                    logger.info(f"Multi-mutation results written to {output_base}")
                except MemoryError:
                    logger.error(
                        "Out of memory during multi-mutation generation. "
                        "Reduce --multi-mutations levels or --multi-top-per-position."
                    )
                    sys.exit(1)
            else:
                write_fasta(args.output, header, sequence,
                            mutations, not args.no_original)
                logger.info(f"Results written to {args.output}")

            print_mutation_summary(mutations)
        else:
            region_analyses = []
            if args.regions:
                for region_str in args.regions:
                    try:
                        start, stop = parse_region(region_str, len(sequence))
                        analysis = analyze_region(
                            sequence, start, stop, selected_rules=args.rules)
                        region_analyses.append(analysis)
                        logger.info(
                            f"Region {start}:{stop}: "
                            f"{len(analysis['merged_clusters'])} clusters found"
                        )
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
        logger.error("Out of memory - reduce --multi-mutations or --multi-top-per-position")
        sys.exit(1)
    except Exception as e:  # pragma: no cover - defensive
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

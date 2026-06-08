"""
FASTA input/output and region-string parsing for AGGRESSOR.

This module is the single source of truth for sequence I/O and region
parsing. Previously these helpers were duplicated across the CLI and the
mutation engine (``parse_region`` / ``_parse_region``), which is a
maintenance hazard: a fix in one copy silently leaves the other stale.
Centralising them here removes that divergence risk.

Functions
---------
read_fasta(path)                         -> (header, sequence)
write_fasta(path, header, seq, muts)     -> None
parse_region("10:50", seq_len)           -> (10, 50)
normalize_regions(["all"], seq_len)      -> ["1:<len>"]
create_output_directory(base, levels)    -> Path
write_multi_mutations_by_category(...)   -> None
"""
from pathlib import Path
from typing import List, Tuple, Dict, Optional

from aggressor.core.config import FASTA_LINE_LENGTH, logger


# =============================================================================
# FASTA READING / WRITING
# =============================================================================

def read_fasta(filepath: str) -> Tuple[str, str]:
    """
    Read a single-sequence FASTA file.

    Only the first record is read (multi-FASTA is not supported); the
    sequence is upper-cased so that downstream residue lookups against the
    canonical 20-letter alphabet are case-insensitive.

    Returns:
        Tuple of (header_line, sequence)

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is malformed or empty.
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
                    break  # single-sequence mode: stop at second record
                header = line
            elif line and not header:
                raise ValueError("FASTA file must start with a header line (>)")
            elif header:
                sequence += line.upper()

        if not header or not sequence:
            raise ValueError("Invalid FASTA file or empty sequence")
        return header, sequence

    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")
    except ValueError:
        raise
    except Exception as e:  # pragma: no cover - defensive
        raise ValueError(f"Error reading FASTA file: {e}")


def _wrap(seq: str) -> List[str]:
    """Wrap a sequence to FASTA_LINE_LENGTH-character lines."""
    return [seq[i:i + FASTA_LINE_LENGTH] for i in range(0, len(seq), FASTA_LINE_LENGTH)]


def write_fasta(
        output_file: str,
        original_header: str,
        original_seq: str,
        mutations: List[Tuple[str, str]],
        include_original: bool = True,
) -> None:
    """Write (description, sequence) records to a FASTA file."""
    protein_name = original_header[1:].strip()
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for chunk in _wrap(original_seq):
                f.write(f"{chunk}\n")
        for description, mutated_seq in mutations:
            f.write(f">{protein_name}_{description}\n")
            for chunk in _wrap(mutated_seq):
                f.write(f"{chunk}\n")


def write_fasta_with_scores(
        output_file: str,
        original_header: str,
        original_seq: str,
        mutations: List[Tuple[str, str, int]],
        include_original: bool = True,
) -> None:
    """Write (description, sequence, score) records; score is not emitted (already in description)."""
    protein_name = original_header[1:].strip()
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for chunk in _wrap(original_seq):
                f.write(f"{chunk}\n")
        for item in mutations:
            description, mutated_seq = item[0], item[1]
            f.write(f">{protein_name}_{description}\n")
            for chunk in _wrap(mutated_seq):
                f.write(f"{chunk}\n")


# =============================================================================
# REGION PARSING
# =============================================================================

def parse_region(region_str: str, seq_length: int) -> Tuple[int, int]:
    """
    Parse a ``start:stop`` region string (1-indexed, inclusive).

    Raises:
        ValueError: If the format is invalid or the region is out of bounds.
    """
    try:
        if ':' not in region_str:
            raise ValueError("Region must be in format start:stop")
        start_str, stop_str = region_str.split(':')
        start = int(start_str.strip())
        stop = int(stop_str.strip())
        if start < 1 or stop > seq_length:
            raise ValueError(f"Region {start}:{stop} out of bounds (1-{seq_length})")
        if start > stop:
            raise ValueError(f"Start position {start} cannot exceed stop position {stop}")
        return start, stop
    except ValueError as e:
        raise ValueError(f"Invalid region format '{region_str}': {e}")


def normalize_regions(
        regions: Optional[List[str]],
        seq_length: int,
) -> Optional[List[str]]:
    """Normalise --regions values; expand the literal 'all' to the full sequence."""
    if not regions:
        return None
    tokens = [r.strip() for r in regions if r is not None and str(r).strip()]
    if not tokens:
        return None
    if any(t.lower() == "all" for t in tokens):
        return [f"1:{seq_length}"]
    return tokens


# =============================================================================
# MULTI-MUTATION OUTPUT STRUCTURE
# =============================================================================

def level_to_text(level: int) -> str:
    """Convert a mutation level to its English name (2 -> 'double')."""
    names = {2: 'double', 3: 'triple', 4: 'quadruple', 5: 'quintuple', 6: 'sextuple'}
    return names.get(level, f'{level}x')


def create_output_directory(
        base_path: str,
        multi_mutation_levels: Optional[List[int]] = None,
) -> Path:
    """Create the organised output directory tree for multi-point mutations."""
    output_dir = Path(base_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    if multi_mutation_levels:
        for level in sorted(multi_mutation_levels):
            level_dir = output_dir / f"{level_to_text(level)}_mutations"
            level_dir.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Created directory: {level_dir}")
    return output_dir


def write_multi_mutations_by_category(
        output_dir: str,
        original_header: str,
        original_seq: str,
        categorized_mutations: Dict[int, Dict[str, List[Dict]]],
        include_original: bool = True,
) -> None:
    """Write categorised multi-point mutations to per-level, per-category FASTA files."""
    output_base = Path(output_dir)
    category_files = [
        ('single_region', 'single_region.fasta'),
        ('multi_region', 'multi_region.fasta'),
        ('all_gatekeeper', 'all_gatekeeper.fasta'),
        ('all_core', 'all_core.fasta'),
        ('mixed', 'mixed.fasta'),
    ]
    for level in sorted(categorized_mutations.keys()):
        level_dir = output_base / f"{level_to_text(level)}_mutations"
        level_dir.mkdir(parents=True, exist_ok=True)
        categories = categorized_mutations[level]
        for category_name, filename in category_files:
            if categories.get(category_name):
                file_path = level_dir / filename
                muts = [
                    (item['description'], item['sequence'], item.get('agg_score', 0))
                    for item in categories[category_name]
                ]
                write_fasta_with_scores(
                    str(file_path), original_header, original_seq, muts, include_original
                )
                logger.info(
                    f"Wrote {len(muts)} {category_name} "
                    f"{level_to_text(level)} mutations to {filename}"
                )

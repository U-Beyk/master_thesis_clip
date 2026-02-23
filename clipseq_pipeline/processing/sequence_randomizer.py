"""
random_genome_sequences.py

Random genome sequence sampler with:

- Exclusion of reference clip sites
- Restriction to gene annotations from GFF3
- Headers kept unchanged
- RAM efficient (indexed genome access)
"""

import random
from collections import defaultdict
from pyfaidx import Fasta


# -------------------------------------------------
# Reference FASTA parsing
# -------------------------------------------------

def parse_header_intervals(fasta_path: str):
    """
    Parse genomic intervals from FASTA headers
    and store headers unchanged.
    """

    intervals = defaultdict(list)
    headers = []

    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            headers.append(header)

            fields = header.split("|")
            field_dict = {}

            for field in fields:
                if ":" in field:
                    key, value = field.split(":", 1)
                    field_dict[key] = value

            clip_range = field_dict["clip_range"]
            chrom = field_dict["chr"]

            start, end = map(int, clip_range.split("-"))
            if end < start:
                start, end = end, start

            intervals[chrom].append((start, end))

    return intervals, headers


# -------------------------------------------------
# GFF3 parsing (gene restriction)
# -------------------------------------------------

def parse_gff3_genes(gff3_path: str):
    """
    Parse GFF3 and return intervals where
    feature type contains 'gene'.
    """

    gene_intervals = defaultdict(list)

    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature_type = parts[2]

            # Include anything containing "gene"
            if "gene" not in feature_type:
                continue

            start = int(parts[3]) - 1  # convert 1-based → 0-based
            end = int(parts[4])        # end stays as-is

            gene_intervals[chrom].append((start, end))

    return gene_intervals


# -------------------------------------------------
# Interval utilities
# -------------------------------------------------

def merge_intervals(intervals):
    merged = {}

    for chrom, ranges in intervals.items():
        if not ranges:
            continue

        ranges = sorted(ranges)
        merged_list = [ranges[0]]

        for start, end in ranges[1:]:
            last_start, last_end = merged_list[-1]

            if start <= last_end:
                merged_list[-1] = (
                    last_start,
                    max(last_end, end),
                )
            else:
                merged_list.append((start, end))

        merged[chrom] = merged_list

    return merged


def overlaps(intervals, chrom, start, end):
    if chrom not in intervals:
        return False

    for s, e in intervals[chrom]:
        if start < e and end > s:
            return True

    return False


def contained_in(intervals, chrom, start, end):
    """
    Check if interval is fully contained
    within at least one allowed interval.
    """

    if chrom not in intervals:
        return False

    for s, e in intervals[chrom]:
        if start >= s and end <= e:
            return True

    return False


# -------------------------------------------------
# Genome handling
# -------------------------------------------------

def load_genome(genome_fasta):
    return Fasta(
        genome_fasta,
        as_raw=True,
        sequence_always_upper=True,
    )


def get_chrom_sizes(genome):
    return {chrom: len(genome[chrom]) for chrom in genome.keys()}


# -------------------------------------------------
# Main function
# -------------------------------------------------

def sample_random_sequences(
    genome_fasta: str,
    reference_fasta: str,
    gff3_file: str,
    seq_length: int,
    output_fasta: str,
    seed: int | None = None,
    max_attempts: int = 1_000_000,
):
    """
    Sample random sequences that:

    - Match number of reference sequences
    - Do NOT overlap reference clip regions
    - Are fully contained within gene annotations
    - Keep original headers unchanged
    """

    if seed is not None:
        random.seed(seed)

    # Parse reference intervals
    exclusion_intervals, headers = parse_header_intervals(reference_fasta)
    exclusion_intervals = merge_intervals(exclusion_intervals)
    n_sequences = len(headers)

    # Parse gene intervals
    gene_intervals = parse_gff3_genes(gff3_file)
    gene_intervals = merge_intervals(gene_intervals)

    genome = load_genome(genome_fasta)
    chrom_sizes = get_chrom_sizes(genome)

    chroms = list(chrom_sizes.keys())
    weights = list(chrom_sizes.values())

    written = 0
    attempts = 0

    with open(output_fasta, "w") as out:

        while written < n_sequences:

            if attempts > max_attempts:
                raise RuntimeError(
                    "Too many sampling attempts. "
                    "Not enough gene space available."
                )

            attempts += 1

            chrom = random.choices(
                chroms, weights=weights, k=1
            )[0]

            max_start = chrom_sizes[chrom] - seq_length
            if max_start <= 0:
                continue

            start = random.randint(0, max_start)
            end = start + seq_length

            # Must be inside gene annotation
            if not contained_in(gene_intervals, chrom, start, end):
                continue

            # Must not overlap exclusion regions
            if overlaps(exclusion_intervals, chrom, start, end):
                continue

            seq = genome[chrom][start:end]

            if "N" in seq:
                continue

            out.write(">" + headers[written] + "\n")
            out.write(seq + "\n")

            written += 1
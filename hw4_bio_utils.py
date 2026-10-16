"""Bioinformatics utility package for sequence filtering and file format conversions."""

"""Main script."""

import os

from helpers.bio_files import (
    extract_top_blast_hits,
    parse_gbk_cds,
    read_multiline_fasta,
)
from helpers.fastq import (
    is_gc_ok,
    is_length_ok,
    is_quality_ok,
    normalize_bounds,
    save_read,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: float | tuple | list = (0, 100),
    length_bounds: float | tuple | list = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """Filter FASTQ reads on the fly by length, GC%, and Phred33 quality."""
    gc_limits = normalize_bounds(gc_bounds)
    length_limits = normalize_bounds(length_bounds)
    min_quality = float(quality_threshold)

    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_fastq)

    with (
        open(input_fastq, "r", encoding="utf-8") as in_f,
        open(output_path, "w", encoding="utf-8") as out_f,
    ):
        while True:
            header = in_f.readline()
            if not header:
                break

            sequence = in_f.readline().strip()
            comment = in_f.readline().strip()
            quality = in_f.readline().strip()
            header = header.strip()

            if not is_length_ok(sequence, length_limits):
                continue
            if not is_gc_ok(sequence, gc_limits):
                continue
            if not is_quality_ok(quality, min_quality):
                continue

            save_read(out_f, header, sequence, comment, quality)


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str | None = None,
) -> None:
    """Convert multiline FASTA file to single-line FASTA."""
    if output_fasta is None:
        base_name, ext = os.path.splitext(input_fasta)
        output_fasta = f"{base_name}_oneline{ext or '.fasta'}"

    records = read_multiline_fasta(input_fasta)

    with open(output_fasta, "w", encoding="utf-8") as out_f:
        out_f.writelines(f"{header}\n{sequence}\n" for header, sequence in records)


def parse_blast_output(input_file: str, output_file: str) -> None:
    """Extract top hit per query from BLAST text file and save sorted."""
    descriptions = extract_top_blast_hits(input_file)
    descriptions.sort()

    with open(output_file, "w", encoding="utf-8") as out_f:
        out_f.writelines(f"{desc}\n" for desc in descriptions)


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: str | list[str] | tuple[str, ...],
    output_fasta: str,
    n_before: int = 1,
    n_after: int = 1,
) -> None:
    """Extract flanking CDS protein translations around target genes into FASTA."""
    if isinstance(genes, str):
        target_genes = {genes}
    else:
        target_genes = set(genes)

    cds_list = parse_gbk_cds(input_gbk)
    total_cds = len(cds_list)
    neighbor_indices = set()

    for i, feature in enumerate(cds_list):
        if feature["gene"] in target_genes:
            start_before = max(0, i - n_before)
            for prev_idx in range(start_before, i):
                neighbor_indices.add(prev_idx)

            end_after = min(total_cds, i + n_after + 1)
            for next_idx in range(i + 1, end_after):
                neighbor_indices.add(next_idx)

    selected_records = [
        cds_list[idx]
        for idx in sorted(neighbor_indices)
        if cds_list[idx]["gene"] not in target_genes
    ]

    with open(output_fasta, "w", encoding="utf-8") as out_f:
        out_f.writelines(
            f">{entry['gene']}\n{entry['translation']}\n" for entry in selected_records
        )

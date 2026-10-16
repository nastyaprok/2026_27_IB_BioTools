"""Parsing utilities for FASTA, BLAST, and GBK files."""


def read_multiline_fasta(file_path: str) -> list[tuple[str, str]]:
    """Read multiline FASTA into list of (header, sequence) pairs."""
    records = []
    current_header = ""
    current_seq_parts = []

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header:
                    records.append((current_header, "".join(current_seq_parts)))
                    current_seq_parts = []
                current_header = line
            else:
                current_seq_parts.append(line)

        if current_header:
            records.append((current_header, "".join(current_seq_parts)))

    return records


def extract_top_blast_hits(file_path: str) -> list[str]:
    """Extract first Description line per query section in BLAST output."""
    best_hits = []
    in_significant_block = False

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if "Sequences producing significant alignments:" in line:
                in_significant_block = True
                continue

            if in_significant_block:
                stripped_line = line.strip()

                if not stripped_line:
                    continue

                if stripped_line.startswith(
                    ("Description", "Scientific", "Max", "Total", "-----")
                ):
                    continue

                if len(line) >= 60 and not line.startswith(" "):
                    description = line[:60].strip()
                else:
                    description = stripped_line.split("  ")[0].strip()

                if description:
                    best_hits.append(description)

                in_significant_block = False

    return best_hits


def parse_gbk_cds(file_path: str) -> list[dict[str, str]]:
    """Parse GBK file into list of CDS dicts with gene and translation."""
    cds_list = []
    is_in_features = False
    current_feature_is_cds = False
    current_gene = None
    current_translation_parts = []
    is_reading_translation = False

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("FEATURES"):
                is_in_features = True
                continue

            if line.startswith(("ORIGIN", "//")):
                if (
                    current_feature_is_cds
                    and current_gene
                    and current_translation_parts
                ):
                    cds_list.append(
                        {
                            "gene": current_gene,
                            "translation": "".join(current_translation_parts),
                        }
                    )
                break

            if not is_in_features:
                continue

            # New feature starts with an indent of 5 spaces
            if line.startswith("     ") and not line.startswith("      "):
                if (
                    current_feature_is_cds
                    and current_gene
                    and current_translation_parts
                ):
                    cds_list.append(
                        {
                            "gene": current_gene,
                            "translation": "".join(current_translation_parts),
                        }
                    )

                feature_type = line.strip().split()[0]
                current_feature_is_cds = feature_type == "CDS"
                current_gene = None
                current_translation_parts = []
                is_reading_translation = False

            elif current_feature_is_cds:
                stripped = line.strip()

                if (
                    stripped.startswith("/gene=")
                    or stripped.startswith("/locus_tag=")
                    and current_gene is None
                ):
                    current_gene = stripped.split("=")[1].strip('"')

                if stripped.startswith("/translation="):
                    is_reading_translation = True
                    stripped = stripped.split("=")[1].lstrip('"')

                if is_reading_translation:
                    current_translation_parts.append(stripped.rstrip('"'))
                    if stripped.endswith('"'):
                        is_reading_translation = False

    return cds_list

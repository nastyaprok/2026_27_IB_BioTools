# Bio.tools Package

A lightweight Python bioinformatics toolkit for FASTQ filtering and parsing sequence annotation files (FASTA, GenBank, BLAST output).

## Requirements
- Python 3.9+
- Only Python Standard Library is used (no third-party dependencies required).

## Installation & Setup
Clone the repository and import modules directly:
```bash
git clone https://github.com/nastyaprok/2026_27_IB_BioTools
cd hw4-utils
```

## Running instructions

To use this bio.tools package, import the required functions into your Python script or interactive interpreter session:

- `from hw4_bio_utils import filter_fastq`
- `from hw4_bio_utils import convert_multiline_fasta_to_oneline`
- `from hw4_bio_utils import parse_blast_output`
- `from hw4_bio_utils import select_genes_from_gbk_to_fasta`

To view docstrings and details about input parameters for any tool, run `help(<function_name>)` in the Python console (e.g., `help(filter_fastq)`).

When running `filter_fastq`, a directory `filtered/` is created automatically in your working path if it does not already exist. All filtered FASTQ files are saved into it under the designated output file name.

The package provides several standalone biological file processing tools:
- `filter_fastq`: Performs on-the-fly streaming filtration of FASTQ reads according to read length, GC percentage, and average Phred33 quality thresholds.
- `convert_multiline_fasta_to_oneline`: Reassembles multiline FASTA sequence records so that each sequence occupies a single continuous line.
- `parse_blast_output`: Scans a BLAST text result file and collects the top hit description for each query, producing an alphabetically sorted single-column text file.
- `select_genes_from_gbk_to_fasta`: Searches a GenBank (`.gbk`) file for designated genes of interest and exports the translated protein sequences of their upstream and downstream neighboring CDS into a FASTA file.
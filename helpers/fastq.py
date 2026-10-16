"""FASTQ helper functions for reading, checking, and saving records."""


def normalize_bounds(bounds: float | tuple | list) -> tuple[float, float]:
    """Convert bounds to (min, max) tuple."""
    if isinstance(bounds, (int, float)):
        return (0.0, float(bounds))
    return (float(bounds[0]), float(bounds[1]))


def is_gc_ok(sequence: str, gc_bounds: tuple[float, float]) -> bool:
    """Check if GC% within given bounds."""
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    gc_percent = (gc_count / len(sequence)) * 100

    min_gc, max_gc = gc_bounds
    return (gc_percent >= min_gc) and (gc_percent <= max_gc)


def is_length_ok(sequence: str, length_bounds: tuple[float, float]) -> bool:
    """Check if read length within given bounds."""
    seq_len = len(sequence)
    min_len, max_len = length_bounds
    return (seq_len >= min_len) and (seq_len <= max_len)


def is_quality_ok(quality_str: str, threshold: float) -> bool:
    """Check if mean Phred33 quality >= threshold."""
    scores = [ord(char) - 33 for char in quality_str]
    mean_quality = sum(scores) / len(scores)
    return mean_quality >= threshold


def save_read(
    file_handle, header: str, sequence: str, comment: str, quality: str
) -> None:
    """Write one FASTQ record to file."""
    file_handle.write(header + "\n")
    file_handle.write(sequence + "\n")
    file_handle.write(comment + "\n")
    file_handle.write(quality + "\n")

# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
"""
Cython-optimized parsing functions for VCF records and string operations.
"""
import re

# Compile regex pattern once at module level
cdef object MOTIF_RANGE_PATTERN = re.compile(r'\((\d+)-(\d+)\)')

def parse_region_fast(str region):
    """
    Fast region parsing: "chr:start-end" -> (chrom, start, end)
    Returns tuple (chrom, start, end) or None if invalid.
    """
    cdef int colon_pos = region.find(':')
    if colon_pos == -1:
        return None
    
    cdef str chrom = region[:colon_pos]
    cdef str positions = region[colon_pos + 1:]
    
    cdef int dash_pos = positions.find('-')
    if dash_pos == -1:
        return None
    
    # Declare variables at function level (required by Cython)
    cdef int start
    cdef int end
    
    try:
        start = int(positions[:dash_pos])
        end = int(positions[dash_pos + 1:])
        return (chrom, start, end)
    except ValueError:
        return None


def format_region_fast(str chrom, int start, int end):
    """
    Fast region formatting: (chrom, start, end) -> "chr:start-end"
    """
    return f"{chrom}:{start}-{end}"


def parse_motif_range_fast(str motif_range):
    """
    Fast motif range parsing using pre-compiled regex.
    Returns list of (start, end) tuples (0-based).
    """
    cdef list matches = MOTIF_RANGE_PATTERN.findall(motif_range)
    cdef list ranges = []
    cdef str start_str, end_str
    cdef int start, end
    
    for start_str, end_str in matches:
        start = int(start_str) - 1
        end = int(end_str) - 1
        ranges.append((start, end))
    
    return ranges


def split_motif_ids_fast(str motif_ids_str, str delimiter="_"):
    """
    Fast splitting of motif IDs string.
    Returns list of motif IDs.
    """
    if not motif_ids_str or motif_ids_str == '.':
        return []
    return motif_ids_str.split(delimiter)


def count_motifs_fast(list motif_ids):
    """
    Fast motif counting using Cython-optimized dictionary operations.
    Returns dict mapping motif_id -> count.
    """
    cdef dict motif_count = {}
    cdef str motif
    cdef int idx
    
    for idx in range(len(motif_ids)):
        motif = motif_ids[idx]
        if motif in motif_count:
            motif_count[motif] += 1
        else:
            motif_count[motif] = 1
    
    return motif_count


def parse_vcf_record_ids_fast(object vcf_file):
    """
    Fast parsing of VCF record IDs and mapping.
    Returns tuple (records_ids dict, records_map dict).
    Optimized for speed with minimal Python overhead.
    """
    cdef dict records_ids = {}
    cdef dict records_map = {}
    cdef int idx = 0
    cdef object rec
    cdef str rec_id, region_str
    
    for rec in vcf_file.fetch():
        rec_id = rec.id
        region_str = f"{rec.chrom}:{rec.pos}-{rec.stop}"
        records_ids[rec_id] = region_str
        records_map[idx] = rec_id
        idx += 1
    
    return records_ids, records_map


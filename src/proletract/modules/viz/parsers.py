"""
Parsing functions for VCF records and motif data processing.
"""
import re
import streamlit as st
import pandas as pd
import pysam

try:
    from ..io.fast_parsing import (
        parse_motif_range_fast as _parse_motif_range_fast,
        count_motifs_fast as _count_motifs_fast,
        split_motif_ids_fast as _split_motif_ids_fast,
        parse_region_fast as _parse_region_fast,
    )
    _USE_CYTHON = True
except (ImportError, ModuleNotFoundError):
    _USE_CYTHON = False
    _parse_motif_range_fast = None
    _count_motifs_fast = None
    _split_motif_ids_fast = None
    _parse_region_fast = None


def parse_record_assembly(vcf, region):
    """
    Extracts a single record from an assembly VCF in the specified region.
    Returns a dictionary with all relevant repeat expansion data.

    Args:
        vcf: A pysam VariantFile or TabixFile object to query or a string (in this case read the vcf file from the path).
        region (str): Region string in the format "chr:start-end".

    Returns:
        dict: All needed fields for downstream processing & visualization.
    """
    if isinstance(vcf, str):
        vcf = pysam.VariantFile(vcf)
        
    chrom, positions = region.split(":")
    try:
        start, end = map(int, positions.split("-"))
    except Exception as e:
        raise ValueError(f"Could not parse genomic region from '{region}' - {e}")
    query_region = f"{chrom}:{start}-{end}"
    try:
        # Try getting the first record in this region
        rec = next(vcf.fetch(region=query_region))
    except StopIteration:
        record = {
            'chr': "",
            'pos': -1,
            'stop': -1,
            'motifs': [],
            'motif_ids_h': [],
            'motif_ids_ref': [],
            'ref_CN': 0,
            'CN_H': 0,
            'spans': [],
            'ref_allele': '',
            'alt_allele': '',
            'gt': '',
            'id': '',
        }
        return record
    # Extract motif ids for the ALT allele
    ids_h = rec.samples[0].get("MI", [])
    if ids_h:
        try:
            if _USE_CYTHON:
                ids_h = _split_motif_ids_fast(ids_h, "_")
            else:
                ids_h = ids_h.split("_")
        except Exception as e:
            st.error(f"Input VCF file is not in the correct format. Please use Reads-based option.")
            st.stop()
    # Extract motif ids for the REF allele
    ids_ref = rec.info.get('MOTIF_IDs_REF', [])
    if ids_ref:
        if _USE_CYTHON:
            ids_ref = _split_motif_ids_fast(ids_ref, "_")
        else:
            ids_ref = ids_ref.split("_")

    # Reference and alternative allele copy numbers
    ref_CN = rec.info.get('CN_ref', 0)
    CN_H = rec.samples[0].get('CN', 0)

    # Get motif names from INFO field
    motif_names = rec.info.get('MOTIFS', [])
    if isinstance(motif_names, tuple):
        motif_names = list(motif_names)
    elif not isinstance(motif_names, list):
        motif_names = [motif_names]

  
    alt_allele = rec.alts[0] if rec.alts and rec.alts[0] != '.' else ''

    # Get the motifs' span
    spans = rec.samples[0].get('SP', [])
    record = {
        'chr': rec.chrom,
        'pos': rec.pos,
        'stop': rec.stop,
        'motifs': motif_names,
        'motif_ids_h': ids_h,
        'motif_ids_ref': ids_ref,
        'ref_CN': ref_CN,
        'CN_H': CN_H,
        'spans': spans,
        'ref_allele': rec.ref,
        'alt_allele': alt_allele,
        'gt': str(rec.samples[0]['GT'][0]),
        'id': rec.id,
    }

    return record


def parse_record(vcf_file, region):
    """
    Parses a single VCF record for a specified region.

    Args:
        vcf_file (str or pysam.VariantFile): VCF file path or pysam file object.
        region (str): Region string in format "chr:start-end".

    Returns:
        dict: Parsed information for the region.
    """
    # Open VCF if a file path is given
    vcf = pysam.VariantFile(vcf_file) if isinstance(vcf_file, str) else vcf_file

    # Fetch the first record in the region
    record_iter = vcf.fetch(region=region)
    rec = next(record_iter, None)
    if rec is None:
        st.warning(f"No records found for region {region}")
        return None

    # Parse motif IDs for both haplotypes
    mi = rec.samples[0]['MI']
    if isinstance(mi, tuple):
        if _USE_CYTHON:
            ids_h1 = _split_motif_ids_fast(mi[0], "_") if mi[0] else []
            ids_h2 = _split_motif_ids_fast(mi[1], "_") if len(mi) > 1 and mi[1] else []
        else:
            ids_h1 = mi[0].split("_") if mi[0] else []
            ids_h2 = mi[1].split("_") if len(mi) > 1 and mi[1] else []
    else:
        if _USE_CYTHON:
            ids_h1 = _split_motif_ids_fast(mi, "_") if mi else []
            ids_h2 = _split_motif_ids_fast(mi, "_") if mi else []
        else:
            ids_h1 = mi.split("_") if mi else []
            ids_h2 = mi.split("_") if mi else []

    # Allele information
    ref_allele = rec.ref
    alt_allele1, alt_allele2 = ".", ""
    if rec.alts:
        alts = list(rec.alts)
        if alts and alts[0] != ".":
            alt_allele1 = alts[0]
        else:
            alt_allele1 = ""
        if len(alts) > 1 and alts[1] != ".":
            alt_allele2 = alts[1]
        elif alts and ids_h1 == ids_h2:
            alt_allele2 = alt_allele1

    # Copy number for both haplotypes
    CNs = list(rec.samples[0]['CN'])
    CN_H1 = str(CNs[0]) if CNs else None
    CN_H2 = str(CNs[1]) if len(CNs) > 1 else None

    # Parse span information
    if 'SP' in rec.samples[0]:
        SP_field = rec.samples[0]['SP']
    else:
        SP_field = "" # in this case tandemtwister genotyped no reads from the sample

    if isinstance(SP_field, tuple):
        spans_h1 = SP_field[0]
        spans_h2 = SP_field[1] if len(SP_field) > 1 else SP_field[0]
        spans = (spans_h1, spans_h2)
    else:
        spans = (SP_field, SP_field)
    # Replace None with empty string
    ref_span = rec.info.get('REF_SPAN', "")
    spans = ["" if x is None else x for x in spans]
    spans = [ref_span] + spans

    # Get motif names from INFO field
    motif_names = rec.info['MOTIFS']
    if isinstance(motif_names, tuple):
        motif_names = list(motif_names)
    elif not isinstance(motif_names, list):
        motif_names = [motif_names]

    # Get genotype from sample field
    gt = rec.samples[0]['GT']
    # Get supporting reads from sample field
    try:
        supporting_reads = rec.samples[0]['DP']
    except:
        st.error(f"Input VCF file is not reads-based VCF files, please use assembly VCF files instead")
        st.stop()
    gt = '/'.join([str(i) for i in gt])
    if isinstance(supporting_reads, tuple):
        supporting_reads_h1 = supporting_reads[0]
        supporting_reads_h2 = supporting_reads[1]
    else:
        supporting_reads_h1 = supporting_reads
        supporting_reads_h2 = supporting_reads
    # Final record dictionary
    record = {
        'chr': rec.chrom,
        'pos': rec.pos,
        'stop': rec.stop,
        'motifs': motif_names,
        'motif_ids_h1': ids_h1,
        'motif_ids_h2': ids_h2,
        'motif_ids_ref': _split_motif_ids_fast(rec.info['MOTIF_IDs_REF'], "_") if _USE_CYTHON else rec.info['MOTIF_IDs_REF'].split("_"),
        'ref_CN': rec.info.get('CN_ref', None),
        'CN_H1': CN_H1,
        'CN_H2': CN_H2,
        'spans': spans,
        'ref_allele': ref_allele,
        'alt_allele1': alt_allele1,
        'alt_allele2': alt_allele2,
        'gt': gt,
        'supported_reads_h1': supporting_reads_h1,
        'supported_reads_h2': supporting_reads_h2,
        'id': rec.id,
    }
    return record


def parse_motif_range(motif_range):
    """Parse motif range string into list of tuples."""
    if _USE_CYTHON:
        return _parse_motif_range_fast(motif_range)
    # Fallback to pure Python
    pattern = re.compile(r'\((\d+)-(\d+)\)')
    matches = pattern.findall(motif_range)
    ranges = [(int(start)-1, int(end)-1) for start, end in matches]
    return ranges


def count_motifs(motif_ids):
    """Count occurrences of each motif ID."""
    if _USE_CYTHON:
        return _count_motifs_fast(motif_ids)
    # Fallback to pure Python
    motif_count = {}
    for idx, motif in enumerate(motif_ids):
        if motif in motif_count:
            motif_count[motif] += 1
        else:
            motif_count[motif] = 1
    return motif_count


def parse_motif_in_region(record, count_motifs_func):
    """
    Parse motif information from a record.
    
    Args:
        record: VCF record dictionary
        count_motifs_func: Function to count motifs (to avoid circular import)
    
    Returns:
        tuple: (motif_names, motif_count_h1, motif_count_h2)
    """
    if record['motif_ids_h1'] == ['.'] and record['motif_ids_h2'] == ['.']:
        return None, None, None

    motif_names = record['motifs']
    motif_count_ref = count_motifs_func(record['motif_ids_ref'])
    found_motifs_ref = list(motif_count_ref.keys())
    found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]
    motif_count_h1 = count_motifs_func(record['motif_ids_h1'])
    found_motifs_h1 = list(motif_count_h1.keys())
    found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
    motif_count_h2 = count_motifs_func(record['motif_ids_h2'])
    found_motifs_h2 = list(motif_count_h2.keys())
    found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
    motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
    return motif_names, motif_count_h1, motif_count_h2


def compute_diploid_genotype_assembly(gt_h1, gt_h2, ids_h1, ids_h2):
    """
    Compute diploid genotype for assembly mode based on haplotypes.
    
    Args:
        gt_h1: Genotype for haplotype 1 (0 or 1, or "0/0", "0/1", etc.)
        gt_h2: Genotype for haplotype 2 (0 or 1, or "0/0", "0/1", etc.)
        ids_h1: Motif IDs for haplotype 1
        ids_h2: Motif IDs for haplotype 2
        
    Returns:
        str: Diploid genotype string (e.g., "0/0", "0/1", "1/1", "1/2")
    """
    # Extract first allele from GT string if it's formatted as "0/1", otherwise use the value directly
    def extract_allele(gt):
        # Handle None or empty string
        if gt is None or gt == '':
            return 0
        if isinstance(gt, str):
            gt = gt.strip()
            if not gt or gt == '.':
                return 0
            if '/' in gt:
                first_allele = gt.split('/')[0].strip()
                if not first_allele or first_allele == '.':
                    return 0
                return int(first_allele)
            else:
                if gt == '.':
                    return 0
                return int(gt)
        else:
            return int(gt) if gt is not None else 0
    
    gt_h1_int = extract_allele(gt_h1)
    gt_h2_int = extract_allele(gt_h2)
    
    # Both are 0 → 0/0
    if gt_h1_int == 0 and gt_h2_int == 0:
        return "0/0"
    
    # One is 0 and one is 1 → 0/1
    elif (gt_h1_int == 0 and gt_h2_int == 1) or (gt_h1_int == 1 and gt_h2_int == 0):
        return "0/1"
    
    # Both are 1 → compare IDs
    elif gt_h1_int == 1 and gt_h2_int == 1:
        # Compare motif IDs - if they're the same, use 1/1, otherwise 1/2 (heterozygous)
        if ids_h1 == ids_h2:
            return "1/1"
        else:
            return "1/2"
    # Fallback - should not happen
    return f"{gt_h1_int}/{gt_h2_int}"


def create_motif_dataframe(sequences, motif_colors, motif_ids, spans_list, motif_names, parse_motif_range_func):
    """
    Create a DataFrame from motif sequences and spans.
    
    Args:
        sequences: List of sequence dictionaries
        motif_colors: Dictionary mapping motif IDs to colors
        motif_ids: List of motif ID lists for each sequence
        spans_list: List of span strings for each sequence
        motif_names: List of motif names
        parse_motif_range_func: Function to parse motif ranges
    
    Returns:
        pd.DataFrame: DataFrame with motif information
    """
    data = []
    interruptions_dict = set()
    for idx, sequence in enumerate(sequences):
        if sequence['sequence'] == "." or sequence['sequence'] == "":
            continue
        sequence_name = sequence['name']
        motif_ids_seq = motif_ids[idx]
        spans = spans_list[idx]
        ranges = parse_motif_range_func(spans)
        sequence_length = len(sequence['sequence'])
        previous_end = 0
        interruptions_dict_sample = {}
        for i, (start, end) in enumerate(ranges):
            motif = motif_ids_seq[i]
            color = motif_colors[int(motif)]

            if start > previous_end:
                data.append({
                    'Sample': sequence_name,
                    'Start': previous_end,
                    'End': start,
                    'Motif': 'Interruption',
                    'Color': '#FF0000',
                    'Sequence': sequence['sequence'][previous_end:start],
                })
                if sequence['sequence'][previous_end:start] in interruptions_dict_sample:
                    interruptions_dict_sample[sequence['sequence'][previous_end:start]] += 1
                else:
                    interruptions_dict_sample[sequence['sequence'][previous_end:start]] = 1
            data.append({
                'Sample': sequence_name,
                'Start': start,
                'End': end + 1,  
                'Motif': motif_names[int(motif)],
                'Color': color,
                'Sequence': sequence['sequence'][start:end+1],
            })

            previous_end = end + 1
        def len_inturruption_is_equal_to_motif_length(motif_names,k):
            for motif in motif_names:
                if len(k) == len(motif):
                    return True
            return False
        
        inturruptions_dict_sample = {k: v for k, v in interruptions_dict_sample.items() if  len_inturruption_is_equal_to_motif_length(motif_names,k) and v > 1}
        for k,v in inturruptions_dict_sample.items():
            interruptions_dict.add(k)
            

        if previous_end < sequence_length:
            data.append({
                'Sample': sequence_name,
                'Start': previous_end,
                'End': sequence_length,
                'Motif': 'Interruption',
                'Color': '#FF0000',
                'Sequence': sequence['sequence'][previous_end:],
            })
    

    return pd.DataFrame(data)


import streamlit as st
import altair as alt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import re 
import pysam
import numpy as np



class Visualization:
    def __init__(self):
        pass


    def get_results_cohort(self,region, files, file_paths):
        samples_results = {}
        chrom, start_end = region.split(":")
        start, end = start_end.split("-")
        start = int(start) - 1
        end = int(end) - 1
        region = f"{chrom}:{start}-{end}"

        for i in range(len(files)):
            sample_name = file_paths[i].split(".")[0]
            record = self.parse_record(files[i], region)
            samples_results[sample_name] = record
        return samples_results
    

        
    def parse_record_assembly(self, vcf, region):
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
        # Parse the chromosome and coordinates, adjust to 0-based for pysam
        chrom, positions = region.split(":")
        try:
            start, end = map(int, positions.split("-"))
        except Exception as e:
            # If region isn't in the correct format, raise a clear error
            raise ValueError(f"Could not parse genomic region from '{region}' - {e}")
        # pysam expects 0-based, half-open intervals
        query_region = f"{chrom}:{start}-{end}"
        try:
            # Try getting the first record in this region. May raise StopIteration.
            rec = next(vcf.fetch(region=query_region))
        except StopIteration:
            # If nothing found, return a 'null record' to indicate missing data
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
            }
            return record
        # Extract motif ids for the ALT allele (usually in sample field)
        ids_h = rec.samples[0].get("MI", [])
        if ids_h:
            ids_h = ids_h.split("_")

        # Extract motif ids for the REF allele (usually in the INFO field)
        ids_ref = rec.info.get('MOTIF_IDs_REF', [])
        if ids_ref:
            ids_ref = ids_ref.split("_")

        # Reference and alternative allele copy numbers
        ref_CN = rec.info.get('CN_ref', 0)
        CN_H = rec.samples[0].get('CN', 0)

        # Get motif names from INFO, and ensure type = list for later use
        motif_names = rec.info.get('MOTIFS', [])
        if isinstance(motif_names, tuple):
            motif_names = list(motif_names)
        elif not isinstance(motif_names, list):
            motif_names = [motif_names]

        # Some VCF encodings use '.' to mean "no ALT"; check for that
        alt_allele = rec.alts[0] if rec.alts and rec.alts[0] != '.' else ''

        # Get the motifs' span, typically from 'SP'
        spans = rec.samples[0].get('SP', [])
        # get the second span
        # Store all relevant information in a single record dict for easy access
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
        }

        return record

    def get_results_hgsvc_pop(self, region, files, file_paths):
        if st.session_state.files == None:
            return None
        samples_results = {}
        for i in range(len(files)):
            sample_name = file_paths[i].split(".")[0]
            record = self.parse_record_assembly(st.session_state.files[i], region)
            samples_results[sample_name] = record
        return samples_results
    

    def visulize_cohort(self):
        col1, middel, col2 = st.columns([1.5,3, 1])  
        if 'cohorts_records_map' in st.session_state:
            region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None, key="region", help="Enter the region in the format: chr:start-end")

            if 'regions_idx' not in st.session_state:
                st.session_state.regions_idx = 0
            with col1:
                if st.button("Previous region"):
                    region = None
                    st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

            with col2:
                if st.button("Next region"):
                    region = None
                    st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len( st.session_state.cohorts_records_map )-1)
            if region and region != st.session_state.get('previous_region', None):
                try:
                    chr_input, start_end_input = region.split(':')
                    start_input, end_input = map(int, start_end_input.split('-'))

                    region = f"{chr_input}:{start_input}-{end_input}"
                    st.session_state.regions_idx = list(st.session_state.cohorts_records_map.values()).index(region)
                except:
                    try:
                        chr_input, start_input, end_input = re.split(r'\s+', region)
                        start_input, end_input = int(start_input), int(end_input)
                        region = f"{chr_input}:{start_input}-{end_input}"
                    except:
                        st.sidebar.info("Invalid region format, showing the first record")
                        region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
                    
                st.session_state.previous_region = region
            else:
                region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
            mode_placeholder = st.empty()

            middel.markdown(f"""
                    <div id="tandem-repeat-region" class="region-container" style="font-size: 25px; margin-bottom: 10px;">
                        <strong>Tandem Repeat Region: {region}</strong>
                    </div>
                """, unsafe_allow_html=True)

            st.markdown("""
                    <style>
                    :root {
                        --region-color-light: black;
                        --region-color-dark: white;
                    }
                    /* Default style for light mode */
                    .region-container {
                        color: var(--region-color-light);
                    }
                    /* Apply different color for dark mode */
                    @media (prefers-color-scheme: dark) {
                        .region-container {
                            color: var(--region-color);
                        }
                    }
                    </style>
                """, unsafe_allow_html=True)

            st.session_state.cohort_results = self.get_results_cohort(region, st.session_state.cohort_files, st.session_state.cohort_file_paths)
            if 'cohort_results' in st.session_state:   
                region = st.session_state.regions_idx
                if (st.session_state.cohort_results == {}):
                    st.warning(f"No motifs found in the region: {region}")
                    st.stop()
                else:
                    self.plot_Cohort_results(st.session_state.cohort_results)
            else:
                st.stop()

    def parse_record_TRGT(self, vcf_file, region):
        """
        Parse a TRGT-formatted VCF record for a given region string (e.g., chr:start-end).
        Handles parsing GT, motif occurrences, alleles, and works out spans for each haplotype.
        Returns a dictionary with all the key info for visualization or downstream analysis.

        Args:
            vcf_file (str): Path to the VCF file.
            region (str): Region string, e.g. "chr1:12345-12399"

        Returns:
            dict: Parsed record info if present in region, otherwise None.
        """



        vcf = pysam.VariantFile(vcf_file)
        record_iter = vcf.fetch(region=region)
        record = next(record_iter, None)

        if record is None:
            st.warning(f"No records found for region {region}")
            return None
        st.write(record)
        # Pull out identifiers and motif info
        ID = record.info['TRID']
        chrom = ID.split('_')[0]
        start = int(ID.split('_')[1])
        end = int(ID.split('_')[2])
        motifs = record.info['MOTIFS']
        sample = record.samples[0]
        GT_tuple = sample['GT']
        GT = '/'.join([str(i) for i in GT_tuple])  

        # fields we'll be filling
        occurrences_h1 = ""
        occurrences_h2 = ""
        motif_ids_hap1 = []
        motif_ids_hap2 = []
        CN_ref = np.nan
        CN_hap1 = np.nan
        CN_hap2 = np.nan
        ref_seq = record.ref
        seq_hap1 = "."
        seq_hap2 = "."

        # Motif "spans" (MS field) -- tells us where the repeats are
        spans = sample.get('MS')
        if spans and spans != ('.',):
            MC_value = sample.get('MC')
            # Support for both tuple and single string for spans
            if isinstance(spans, tuple) and len(spans) > 0:
                spans_string_h1 = spans[0] if spans[0] else ""
                spans_string_h2 = spans[1] if len(spans) > 1 else spans[0]
                MC_value_h1 = MC_value[0] if MC_value[0] else ""
                MC_value_h2 = MC_value[1] if len(MC_value) > 1 else MC_value[0]
                MC_value_h1 = list(map(int, MC_value_h1.split('_')))
                MC_value_h2 = list(map(int, MC_value_h2.split('_')))
            else:
                spans_string_h1 = spans
                spans_string_h2 = spans
                MC_value_h1 = MC_value if MC_value else ""
                MC_value_h2 = MC_value if MC_value else ""

            motif_indices_h1 = []
            motif_indices_h2 = []
            intervals_h1 = []
            intervals_h2 = []
            # Split and parse the motif indices and intervals for both haplotypes.
            def parse_spans(spans_string, motif_indices, intervals):
                if spans_string:
                    spans_list = spans_string.split('_')
                    a = 0
                    for span in spans_list:
                        match = re.match(r'(\d+)\((\d+)-(\d+)\)', span)
                        if match:
                            motif_indices.append(str(match.group(1)))
                            # Store tuple of ints for easier downstream adjustment
                            if a==0:
                                intervals.append((int(match.group(2))+1, int(match.group(3))+1))
                            else:
                                intervals.append((int(match.group(2))+2, int(match.group(3))+1))
                            a += 1
            parse_spans(spans_string_h1, motif_indices_h1, intervals_h1)
            parse_spans(spans_string_h2, motif_indices_h2, intervals_h2)
            # make the intervals in this format again (start-end    )
            intervals_h1 = [f"({start}-{end})" for start, end in intervals_h1]
            intervals_h2 = [f"({start}-{end})" for start, end in intervals_h2]
 
            spans_h1 = ''.join(intervals_h1)
            spans_h2 = ''.join(intervals_h2)
   
            # Get occurrence numbers from spans strings
            # This function must be implemented elsewhere in the class
            # occurrences_h1 = self.convert_trgt_spans_to_occurrences(spans_string_h1, motifs)
            # occurrences_h2 = self.convert_trgt_spans_to_occurrences(spans_string_h2, motifs)

            # # For each haplotype, the occurrences string is like "3_2" (for motif counts)
            # motif_ids_hap1 = list(map(int, occurrences_h1.split('_'))) if occurrences_h1 else []
            # motif_ids_hap2 = list(map(int, occurrences_h2.split('_'))) if occurrences_h2 else []
            # get the MC format
            CN_hap1 = sum(MC_value_h1)
            CN_hap2 = sum(MC_value_h2)

            # # Convert motif counts to string format for output
            # motif_ids_hap1_str = [str(i) for i in MC_value_h1]
            # motif_ids_hap2_str = [str(i) for i in MC_value_h2]

      

            # Figure out alternate allele sequences
            alt_allele1 = "."
            alt_allele2 = "."

            if record.alts:
                # Usually one or two alternates; handle accordingly
                if len(record.alts) > 0:
                    alt_allele2 = record.alts[0]
                    if len(record.alts) > 1:
                        alt_allele2 = record.alts[1]
                    else:
                        # Only one alt; what haplotype does it belong to? 
                        if GT == '0/1':
                            alt_allele1 = ref_seq
                        else:
                            # Both haplotypes are alt, GT==1/1
                            alt_allele1 = alt_allele2
                else:
                    # No valid alt, just reference
                    alt_allele1 = ref_seq
                    alt_allele2 = ref_seq
            else:
                if GT == '0/0':
                    alt_allele1 = ref_seq
                    alt_allele2 = ref_seq
                
            # Set output sequences for haplotypes if available
            if alt_allele1 != ".":
                seq_hap1 = alt_allele1
            if alt_allele2 != ".":
                seq_hap2 = alt_allele2


            motif_ids_hap1_str = motif_indices_h1
            motif_ids_hap2_str = motif_indices_h2
            # # Split the overall spans for visualization (especially for stacked motif visualization)
            # spans_h1_intervals = split_spans_by_occurrences(spans_string_h1, motif_ids_hap1, motifs)
            # spans_h2_intervals = split_spans_by_occurrences(spans_string_h2, motif_ids_hap2, motifs)
            # spans_list = [spans_h1_intervals, spans_h2_intervals]
            spans_list = [spans_h1, spans_h2]
        else:
            motif_ids_hap1_str = []
            motif_ids_hap2_str = []
            spans_list = []

        output_record = {
            'chr': chrom,
            'pos': start,
            'stop': end,
            'motifs': list(motifs),
            'motif_ids_h1': motif_ids_hap1_str,
            'motif_ids_h2': motif_ids_hap2_str,
            'motif_ids_ref': [],
            'ref_CN': CN_ref,
            'CN_H1': CN_hap1,
            'CN_H2': CN_hap2,
            'spans': spans_list,
            'ref_allele': ref_seq,
            'alt_allele1': seq_hap1,
            'alt_allele2': seq_hap2
        }
        return output_record

    def parse_record(self, vcf_file, region):
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

        # Parse motif IDs for both haplotypes (MI field)
        mi = rec.samples[0]['MI']
        if isinstance(mi, tuple):
            ids_h1 = mi[0].split("_") if mi[0] else []
            ids_h2 = mi[1].split("_") if len(mi) > 1 and mi[1] else []
        else:
            ids_h1 = mi.split("_") if mi else []
            ids_h2 = mi.split("_") if mi else []

        # Allele information
        ref_allele = rec.ref
        alt_allele1, alt_allele2 = ".", ""
        if rec.alts:
            alts = list(rec.alts)
            # Assign first allele
            if alts and alts[0] != ".":
                alt_allele1 = alts[0]
            else:
                alt_allele1 = ""
            # Assign second allele
            if len(alts) > 1 and alts[1] != ".":
                alt_allele2 = alts[1]
            elif alts and ids_h1 == ids_h2:
                alt_allele2 = alt_allele1

        # Copy number for both alleles
        CNs = list(rec.samples[0]['CN'])
        CN_H1 = str(CNs[0]) if CNs else None
        CN_H2 = str(CNs[1]) if len(CNs) > 1 else None

        # Parse span information (SP field), fallback to empty if not present
        if 'SP' in rec.samples[0]:
            SP_field = rec.samples[0]['SP']
        else:
            SP_field = "" # in this case tandemtwister genotyped no reads from the sample (both CNs are 0, DELETED)

        if isinstance(SP_field, tuple):
            spans_h1 = SP_field[0]
            spans_h2 = SP_field[1] if len(SP_field) > 1 else SP_field[0]
            spans = (spans_h1, spans_h2)
        else:
            spans = (SP_field, SP_field)
        # Replace None with empty string, prepend reference span
        ref_span = rec.info.get('REF_SPAN', "")
        spans = ["" if x is None else x for x in spans]
        spans = [ref_span] + spans

        # Parse motif names
        motif_names = rec.info['MOTIFS']
        if isinstance(motif_names, tuple):
            motif_names = list(motif_names)
        elif not isinstance(motif_names, list):
            motif_names = [motif_names]

        # Final record dictionary
        record = {
            'chr': rec.chrom,
            'pos': rec.pos,
            'stop': rec.stop,
            'motifs': motif_names,
            'motif_ids_h1': ids_h1,
            'motif_ids_h2': ids_h2,
            'motif_ids_ref': rec.info['MOTIF_IDs_REF'].split("_"),
            'ref_CN': rec.info.get('CN_ref', None),
            'CN_H1': CN_H1,
            'CN_H2': CN_H2,
            'spans': spans,
            'ref_allele': ref_allele,
            'alt_allele1': alt_allele1,
            'alt_allele2': alt_allele2
        }
        # For debugging/inspection purposes
        return record

    def compare_different_technologies(self):
        if 'regions_idx' not in st.session_state:
            st.session_state.regions_idx = 0
        st.sidebar.markdown("### Select Region to Visualize")
        region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None, key="region", help="Enter the region in the format: chr:start-end")
        col1, middel, col2 = st.columns([1.5,3, 1])  
        REF, CN1_col, CN2_col = st.columns([1, 1, 1])
        with col1:
            if st.button("Previous region"):
                region = None
                st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)
        with col2:
            if st.button("Next region"):
                region = None
                st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.assembly_vcf_records_map_h1) - 1)
        if region and region != st.session_state.get('previous_region', None):
            try:
                chr_input, start_end_input = region.split(':')
                start_input, end_input = map(int, start_end_input.split('-'))
                input_region = f"{chr_input}:{start_input}-{end_input}"
                if input_region not in st.session_state.assembly_vcf_records_h1:
                    st.warning(f"No records found for the region: {input_region}")
                record_key = st.session_state.assembly_vcf_records_h1[input_region]
                st.session_state.regions_idx = list(st.session_state.assembly_vcf_records_map_h1.values()).index(input_region)

            except:
                try:
                    chr_input, start_input, end_input = re.split(r'\s+', region)
                    start_input, end_input = int(start_input), int(end_input)
                    input_region = f"{chr_input}:{start_input}-{end_input}"
                    record_key = st.session_state.assembly_vcf_records_h1[input_region]
                    st.session_state.regions_idx = list(st.session_state.assembly_vcf_records_map_h1.values()).index(input_region)
                except:
                    st.warning(f"Invalid region format: {region}")
                    st.stop()
        else:
            try:
                record_key = st.session_state.assembly_vcf_records_h1[st.session_state.assembly_vcf_records_map_h1[st.session_state.regions_idx]]
            except:
                try:
                    record_key = st.session_state.assembly_vcf_records_h1[st.session_state.assembly_vcf_records_map_h1[st.session_state.regions_idx]]   
                except:
                    st.warning(f"No records found for the region in the assembly vcfs")
                    st.stop()
                    
        st.session_state.previous_region = region
        tandemtwister_record = self.parse_record(st.session_state.vcf_file_tandemtwister, record_key)
        assembly_record_h1 = self.parse_record_assembly(st.session_state.assembly_vcf_h1, record_key)
        assembly_record_h2 = self.parse_record_assembly(st.session_state.assembly_vcf_h2, record_key)
        trgt_record = self.parse_record_TRGT(st.session_state.vcf_file_trgt, record_key)



        # /confidential/home01/Calraei/tandemrepeats/trgt/results/HG002_trgt.vcf.gz
        # /confidential/home01/Calraei/tandemrepeats/tandemtwister/results/HG002_CCS.vcf.gz
        # /confidential/home01/Calraei/tandemrepeats/tandemtwister/assembly_results/NA24385_h1.vcf.gz
        if trgt_record is None or tandemtwister_record is None or assembly_record_h1 is None or assembly_record_h2 is None:
            st.warning(f"No records found for the region: {record_key}")
            st.stop()
        st.markdown(f"""
                    <div id="tandem-repeat-region" class="region-container" style="font-size: 25px; margin-bottom: 10px;">
                        <strong>Tandem Repeat Region: {record_key}</strong>
                    </div>
                """, unsafe_allow_html=True)
        
        st.markdown("""
                <style>
                :root {
                    --region-color-light: black;
                    --region-color-dark: white;
                }
                /* Default style for light mode */
                .region-container {
                    color: var(--region-color-light);
                }
                /* Apply different color for dark mode */
                @media (prefers-color-scheme: dark) {
                    .region-container {
                        color: var(--region-color);
                    }
                }
                </style>
            """, unsafe_allow_html=True)

        self.compare_TR_results(trgt_record, tandemtwister_record, assembly_record_h1, assembly_record_h2)
    def merge_spans_based_on_motif_ids(self, spans, motif_ids):
        """
        Given a list of spans (as strings) and motif IDs, merges consecutive spans
        with the same motif ID and returns a list of merged span strings and the corresponding merged ids.
        
        Args:
            spans (list): List of strings like '(0-8)', '(8-12)', ...
            motif_ids (list): List of motif IDs, e.g., [0,0,0,2,1]
        Returns:
            Tuple:
              - List of merged spans as strings, e.g. ['(0-16)', '(16-19)', '(19-23)']
              - List of the merged motif IDs, e.g. [0, 2, 1]
        """
        import re

        def parse_span(span_str):
            # Expect input like '(0-8)' or similar
            s = str(span_str).strip()
            if not (s.startswith('(') and s.endswith(')')):
                return None, None
            s = s[1:-1]  # Remove parentheses
            parts = s.split('-')
            if len(parts) != 2:
                return None, None
            try:
                return int(parts[0]), int(parts[1])
            except ValueError:
                return None, None


        # the span is a str  (1-3)(4-6)(7-9)(10-12)(13-15)(16-18)(19-21)(22-26)(27-31)(32-35)(36-39)(40-43)(44-47)(48-50)(51-54)
        # split the spans into a list of spans
        # add _ after each span
        tmp_span_str = ""
        for i in range(len(spans)):
            if spans[i] == ")":
                tmp_span_str += spans[i] + "_"
            else:
                tmp_span_str += spans[i]
        spans = tmp_span_str.split("_")
        spans = [s for s in spans if s != ""]

        if not spans or not motif_ids or len(spans) != len(motif_ids):
            return ".", ["."]
        merged_spans = []
        merged_ids = []
        curr_start, curr_end = None, None
        curr_id = None

        prev_end = None
        prev_seq_int = None

        for idx, (span, mid) in enumerate(zip(spans, motif_ids)):
            s, e = parse_span(span)
            if s is None or e is None:
                continue
            # At first, or after a merge, set tracking variables
            if curr_id is None:
                curr_id = mid
                curr_start = s
                curr_end = e
                prev_end = e
                prev_seq_int = None if idx == 0 else None
            else:
                # Check if motifs match and intervals are contiguous (no interruption/gap)
                # That is, if this span starts exactly after the last one ended
                if mid == curr_id and s == prev_end + 1:
                    # Extend current span
                    curr_end = e
                    prev_end = e
                else:
                    # End current, start new
                    merged_spans.append(f"({curr_start}-{curr_end})")
                    merged_ids.append(curr_id)
                    curr_id = mid
                    curr_start = s
                    curr_end = e
                    prev_end = e
        # Don't forget the last run
        if curr_start is not None and curr_end is not None:
            merged_spans.append(f"({curr_start}-{curr_end})")
            merged_ids.append(curr_id)
        merged_spans = "".join(merged_spans)
        return merged_spans, merged_ids
    def visulize_assembly_h1_h2(self, assembly_record_h1, assembly_record_h2):
     
        # Visualize the dynamic motifs for assembly haplotype 1 and 2
        assembly_h1_motifs_ids = assembly_record_h1['motif_ids_h'] if assembly_record_h1['motif_ids_h'] != ['.'] else []
        assembly_h2_motifs_ids = assembly_record_h2['motif_ids_h'] if assembly_record_h2['motif_ids_h'] != ['.'] else []

    
        assembly_h1_spans = assembly_record_h1['spans'] if len(assembly_record_h1['spans']) > 0 else "."
        assembly_h2_spans = assembly_record_h2['spans'] if len(assembly_record_h2['spans']) > 0 else "."
 
        #  sequence_name, sequence, motif_ids, spans, motif_colors, motif_names):
        assembly_h1_motifs_names = assembly_record_h1['motifs'] if len(assembly_record_h1['motifs']) > 0 else "."
        assembly_h2_motifs_names = assembly_record_h2['motifs'] if len(assembly_record_h2['motifs']) > 0 else "."

        assembly_h1_alt_allele = assembly_record_h1['alt_allele'] if len(assembly_record_h1['alt_allele']) > 0 else "."
        assembly_h2_alt_allele = assembly_record_h2['alt_allele'] if len(assembly_record_h2['alt_allele']) > 0 else "."

        assembly_h1_spans, assembly_h1_motifs_ids = self.merge_spans_based_on_motif_ids(assembly_h1_spans, assembly_h1_motifs_ids)
        assembly_h2_spans, assembly_h2_motifs_ids = self.merge_spans_based_on_motif_ids(assembly_h2_spans, assembly_h2_motifs_ids)
 
        assembly_h1_motifs_colors = self.get_color_palette(len(assembly_h1_motifs_names))
        assembly_h1_motifs_colors = {idx: color for idx, color in enumerate(assembly_h1_motifs_colors)}
        assembly_h2_motifs_colors = self.get_color_palette(len(assembly_h2_motifs_names))
        assembly_h2_motifs_colors = {idx: color for idx, color in enumerate(assembly_h2_motifs_colors)}
        if assembly_h1_motifs_names == ['.'] and assembly_h2_motifs_names == ['.']:
            st.warning(f"No motifs found in the region: {assembly_record_h1['chr']}:{assembly_record_h1['pos']}-{assembly_record_h1['stop']}")
            st.stop()
        else:
            if assembly_h1_motifs_names == ['.']:
                assembly_h1_motifs_names = assembly_h2_motifs_names
                assembly_h1_motifs_colors = assembly_h2_motifs_colors

            else:
                assembly_h2_motifs_names = assembly_h1_motifs_names
                assembly_h2_motifs_colors = assembly_h1_motifs_colors

        
        self.display_dynamic_sequence_with_highlighted_motifs(
            "AS H1",
            assembly_h1_alt_allele,        
            assembly_h1_motifs_ids,      
            assembly_h1_spans,                        
            assembly_h1_motifs_colors,   
            assembly_h1_motifs_names    
        )
        self.display_dynamic_sequence_with_highlighted_motifs(
            "AS H2",
            assembly_h2_alt_allele,        
            assembly_h2_motifs_ids,     
            assembly_h2_spans,       
            assembly_h2_motifs_colors,   
            assembly_h2_motifs_names   
        )


        
    def visulize_tandemtwister_h1_h2(self, tandemtwister_record, tandemtwister_motifs_colors):

        tandemtwister_motifs_names = tandemtwister_record['motifs']
        tandemtwister_motifs_ids_h1 = tandemtwister_record['motif_ids_h1'] if tandemtwister_record['motif_ids_h1'] != ['.'] else []
        tandemtwister_motifs_ids_h2 = tandemtwister_record['motif_ids_h2'] if tandemtwister_record['motif_ids_h2'] != ['.'] else []
        tandemtwister_spans_h1 = tandemtwister_record['spans'][1]  if len(tandemtwister_record['spans']) > 1 else "."
        tandemtwister_spans_h2 = tandemtwister_record['spans'][2] if len(tandemtwister_record['spans']) > 2 else "."
        tandemtwister_alt_allele_h1 = tandemtwister_record['alt_allele1']
        tandemtwister_alt_allele_h2 = tandemtwister_record['alt_allele2']
        
        tandemtwister_spans_h1, tandemtwister_motifs_ids_h1 = self.merge_spans_based_on_motif_ids(tandemtwister_spans_h1, tandemtwister_motifs_ids_h1)
        tandemtwister_spans_h2, tandemtwister_motifs_ids_h2 = self.merge_spans_based_on_motif_ids(tandemtwister_spans_h2, tandemtwister_motifs_ids_h2)

        self.display_dynamic_sequence_with_highlighted_motifs(
                "TT H1",
                tandemtwister_alt_allele_h1,        
                tandemtwister_motifs_ids_h1,     
                tandemtwister_spans_h1,       
                tandemtwister_motifs_colors,   
                tandemtwister_motifs_names   
            )
        self.display_dynamic_sequence_with_highlighted_motifs(
            "TT H2",
            tandemtwister_alt_allele_h2,        
            tandemtwister_motifs_ids_h2,     
            tandemtwister_spans_h2,       
            tandemtwister_motifs_colors,   
            tandemtwister_motifs_names   
        )
    def visulize_trgt_results(self, trgt_record, trgt_motifs_colors):
        trgt_motifs_names = trgt_record['motifs']
        trgt_motifs_ids_h1 = trgt_record['motif_ids_h1'] 
        trgt_motifs_ids_h2 = trgt_record['motif_ids_h2']
        trgt_spans_h1 = trgt_record['spans'][0] if len(trgt_record['spans']) > 0 else "."
        trgt_spans_h2 = trgt_record['spans'][1] if len(trgt_record['spans']) > 1 else "."
        trgt_alt_allele_h1 = trgt_record['alt_allele1'] if len(trgt_record['alt_allele1']) > 0 else "."
        trgt_alt_allele_h2 = trgt_record['alt_allele2'] if len(trgt_record['alt_allele2']) > 0 else "."
  
        self.display_dynamic_sequence_with_highlighted_motifs(
            "TG H1",
            trgt_alt_allele_h1,
            trgt_motifs_ids_h1,
            trgt_spans_h1,
            trgt_motifs_colors,
            trgt_motifs_names
        )
        self.display_dynamic_sequence_with_highlighted_motifs(
            "TG H2",
            trgt_alt_allele_h2,        
            trgt_motifs_ids_h2,     
            trgt_spans_h2,       
            trgt_motifs_colors,   
            trgt_motifs_names   
        )
    def compare_TR_results(self, trgt_record, tandemtwister_record, assembly_record_h1, assembly_record_h2):
        # print the motif ids and the spans of the TRGT, TandemTwister, and Assembly records
        assembly_motifs_names = assembly_record_h1['motifs']
        trgt_motifs_names = trgt_record['motifs']
        tandemtwister_motifs_names = tandemtwister_record['motifs']
        # make a color map and assign it to the motifs but only for the assembly motifs
        assembly_motifs_colors = self.get_color_palette(len(assembly_motifs_names))
        assembly_motifs_colors = {assembly_motif: color for assembly_motif, color in zip(assembly_motifs_names, assembly_motifs_colors)}
        # based on the assembly motifs colors, assign the colors to the trgt and tandemtwister motifs
        trgt_motifs_colors = {idx: assembly_motifs_colors[trgt_motif] for idx, trgt_motif in enumerate(trgt_motifs_names)}
        tandemtwister_motifs_colors = {idx: assembly_motifs_colors[tandemtwister_motif] for idx, tandemtwister_motif in enumerate(tandemtwister_motifs_names)}
        assembly_motifs_colors = {idx: assembly_motifs_colors[assembly_motif] for idx, assembly_motif in enumerate(assembly_motifs_names)}
  
        # visulize the assembly h1 and h2 with coloring the spans with the motif ids
        col1 , col_middle, col2 = st.columns([3,0.2,0.5])
        with col1:
            self.visulize_assembly_h1_h2(assembly_record_h1, assembly_record_h2)
            self.visulize_tandemtwister_h1_h2(tandemtwister_record, tandemtwister_motifs_colors)
            self.visulize_trgt_results(trgt_record, trgt_motifs_colors)
        with col2:
            self.display_motif_legend(tandemtwister_motifs_names, tandemtwister_motifs_colors,col2)

        # Visualize a table showing the CNs (copy numbers) of the TRGT, TandemTwister, and Assembly records

        def get_allele_CNs(record):
            # Utility to get allele CN info from a record, handling missing fields gracefully
            
            cn1 = record.get('CN_H1', None)
            cn2 = record.get('CN_H2', None)
            # If CN fields missing, try to infer from spans (by counting '-' pairs as a proxy)
            if cn1 is None and 'spans' in record and len(record['spans']) > 1:
                cn1 = str(record['spans'][0]).count('-')
            if cn2 is None and 'spans' in record and len(record['spans']) > 1:
                cn2 = str(record['spans'][1]).count('-')
            return cn1, cn2
        def get_assembly_CNs(record):
            cn = record.get('CN_H', None)
            return cn
        trgt_cn1, trgt_cn2 = get_allele_CNs(trgt_record)
        tandemtwister_cn1, tandemtwister_cn2 = get_allele_CNs(tandemtwister_record)
        assembly_h1_cn = get_assembly_CNs(assembly_record_h1)
        assembly_h2_cn = get_assembly_CNs(assembly_record_h2)

        data = {
            "Caller": ["TRGT", "TandemTwister", "Assembly"],
            "Allele 1 CN": [trgt_cn1, tandemtwister_cn1, assembly_h1_cn],
            "Allele 2 CN": [trgt_cn2, tandemtwister_cn2, assembly_h2_cn]
        }
        import pandas as pd
        cn_table = pd.DataFrame(data)
        st.markdown("### CNs (Copy Numbers) from Each Tool/Allele")
        st.dataframe(cn_table, hide_index=True)
        
        
        
      
    def visulize_region(self):
        if 'regions_idx' not in st.session_state:
            st.session_state.regions_idx = 0
        st.sidebar.markdown("### Select Region to Visualize")
        region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None, key="region", help="Enter the region in the format: chr:start-end")
                
        display_option = st.sidebar.radio("Select Display Type", 
                                            ("Sequence with Highlighted Motifs", "Bars"))

        col1, middel, col2 = st.columns([1.5,3, 1])  
        REF, CN1_col, CN2_col = st.columns([1, 1, 1])
        with col1:
            if st.button("Previous region"):
                region = None
                st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

        with col2:
            if st.button("Next region"):
                region = None
                st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.records_map) - 1)
        if region and region != st.session_state.get('previous_region', None):
            try:
                chr_input, start_end_input = region.split(':')
                start_input, end_input = map(int, start_end_input.split('-'))
           

                input_region = f"{chr_input}:{start_input}-{end_input}"
                record_key = st.session_state.records[input_region]
                st.session_state.regions_idx = list(st.session_state.records_map.values()).index(input_region)
                    

            except:
                try:
                    chr_input, start_input, end_input = re.split(r'\s+', region)
                    start_input, end_input = int(start_input), int(end_input)
                    input_region = f"{chr_input}:{start_input}-{end_input}"
                    record_key = st.session_state.records[input_region]
                    st.session_state.regions_idx = list(st.session_state.records_map.values()).index(input_region)
                except:
                    st.sidebar.info("Invalid region format, showing the first record")
                    record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]
        else:
            record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]
        st.session_state.previous_region = region
        record = self.parse_record(st.session_state.vcf_file_path, record_key)
        if record is None:
            st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
            st.stop()
        hgsvc_records = self.get_results_hgsvc_pop(record_key, st.session_state.files ,st.session_state.file_paths)
        if len(record["motif_ids_h1"]) == 0 and len(record["motif_ids_h2"]) == 0:
            st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
            st.stop()
                    
        middel.markdown(f"""
                        <div id="tandem-repeat-region" class="region-container" style="font-size: 25px; margin-bottom: 10px;">
                            <strong>Tandem Repeat Region: {st.session_state.records_map[st.session_state.regions_idx]}</strong>
                        </div>
                    """, unsafe_allow_html=True)

        st.markdown("""
                        <style>
                        :root {
                            --region-color-light: black;
                            --region-color-dark: white;
                        }
                        /* Default style for light mode */
                        .region-container {
                            color: var(--region-color-light);
                        }
                        /* Apply different color for dark mode */
                        @media (prefers-color-scheme: dark) {
                            .region-container {
                                color: var(--region-color);
                            }
                        }
                        </style>
                    """, unsafe_allow_html=True)
        REF.markdown(f"""
                    <div style="font-size: 20px; color: #4CAF50; margin-bottom: 10px;">
                        <strong>Reference Copy Number:</strong> {record['ref_CN']}
                    </div>
                """, unsafe_allow_html=True)

        left_column, right_column = st.columns([4, 1])
        motif_colors = self.get_color_palette(len(record['motifs']))
        motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
            
        col1,col2 = st.sidebar.columns([1,1])



        if display_option == "Sequence with Highlighted Motifs":
            self.visulize_TR_with_dynamic_sequence(record,hgsvc_records, left_column, right_column,motif_colors,CN1_col,CN2_col, st.session_state.get('show_comparison', False))

        elif display_option == "Bars":
            self.display_motifs_with_bars(record, left_column, right_column,motif_colors,CN1_col,CN2_col, st.session_state.get('show_comparison', False),hgsvc_records)


    def parse_motif_range(self,motif_range):
        pattern = re.compile(r'\((\d+)-(\d+)\)')
        matches = pattern.findall(motif_range)
        ranges = [(int(start)-1, int(end)-1) for start, end in matches]
        return ranges
    
    def convert_trgt_spans_to_occurrences(self, spans_string, motifs):
        """
        Convert TRGT spans string to motif occurrences based on motif sizes and spans.
        
        Args:
            spans_string (str): TRGT spans format like '0(3-9)_2(15-21)_3(22-29)_4(29-35)'
            motifs (tuple/list): List of motif sequences to determine motif lengths
            
        Returns:
            str: Motif occurrences like '0_2_3_4' or '0_0_2_3_4' if first motif spans longer
        """
        if not spans_string or spans_string == "":
            return ""
            
        # Parse the spans string to extract motif index and position ranges
        # Format: '0(3-9)_2(15-21)_3(22-29)_4(29-35)'
        pattern = re.compile(r'(\d+)\((\d+)-(\d+)\)')
        matches = pattern.findall(spans_string)
        
        occurrences = []
        for motif_idx_str, start_str, end_str in matches:
            motif_idx = int(motif_idx_str)
            start_pos = int(start_str)
            end_pos = int(end_str)
           
            # Get motif length from the motifs tuple/list
            if motif_idx < len(motifs):
                motif_length = len(motifs[motif_idx])
            else:
                # Fallback if motif index is out of range
                motif_length = end_pos - start_pos + 1
        
            # Calculate how many motif occurrences fit in this span
            span_length = end_pos - start_pos + 1
            num_occurrences = span_length // motif_length
            # Add the motif index the appropriate number of times
            for _ in range(num_occurrences):
                occurrences.append(str(motif_idx))
        
        return '_'.join(occurrences)   
    
    def parse_motif_in_region(self,record):

        if record['motif_ids_h1'] == ['.'] and record['motif_ids_h2'] == ['.']:
            return None, None, None

        motif_names = record['motifs']
        motif_count_ref = self.count_motifs(record['motif_ids_ref'])
        found_motifs_ref = list(motif_count_ref.keys())
        found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]
        motif_count_h1 = self.count_motifs(record['motif_ids_h1'])
        found_motifs_h1 = list(motif_count_h1.keys())
        found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
        motif_count_h2 = self.count_motifs(record['motif_ids_h2'])
        found_motifs_h2 = list(motif_count_h2.keys())
        found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
        motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
        motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
        return motif_names,motif_count_h1,motif_count_h2

    # def display_dynamic_sequence_with_highlighted_motifs(self, sequence_name, sequence, motif_ids, spans, motif_colors, motif_names):
    #     # Handle the situations where motifs are not available or sequence is just one base

    #     if (motif_ids == ["."]) or (isinstance(sequence, str) and len(sequence) <= 1):
    #         if sequence_name == "Ref":
    #             sequence_name += "seq"
    #         st.markdown(f"""
    #             <style>
    #                 .sequence-container {{
    #                     font-family: 'SF Mono', 'Monaco', 'Consolas', 'Roboto Mono', monospace;
    #                     margin-bottom: 20px;
    #                     border: 1px solid #e1e5e9;
    #                     border-radius: 12px;
    #                     overflow: hidden;
    #                     box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    #                     background: white;
    #                 }}
    #                 .sequence-header {{
    #                     background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    #                     color: white;
    #                     padding: 12px 16px;
    #                     font-weight: 600;
    #                     font-size: 14px;
    #                     display: flex;
    #                     justify-content: space-between;
    #                     align-items: center;
    #                 }}
    #                 .sequence-length {{
    #                     font-size: 12px;
    #                     opacity: 0.9;
    #                     font-weight: normal;
    #                 }}
    #                 .no-motif-msg {{
    #                     padding: 16px;
    #                     background: #fff7f7;
    #                     color: #e53935;
    #                     font-weight: 500;
    #                     border-radius: 8px;
    #                     border: 1px solid #f1c0c0;
    #                     margin-top: 10px;
    #                     margin-bottom: 10px;
    #                     font-size: 15px;
    #                     text-align: left;
    #                 }}
    #                 .sequence-content {{
    #                     padding: 16px;
    #                     max-width: 100%;
    #                     white-space: nowrap;
    #                     overflow-x: auto;
    #                     overflow-y: hidden;
    #                     background: #f8f9fa;
    #                     border-top: 1px solid #e1e5e9;
    #                     line-height: 1.6;
    #                     scrollbar-width: thin;  /* Firefox */
    #                     scrollbar-color: #b5b5b5 #f8f9fa; /* Firefox */
    #                 }}
    #             </style>
    #             <div class="sequence-container">
    #                 <div class="sequence-header">
    #                     <span>{sequence_name}</span>
    #                     <span class="sequence-length">Length: {len(sequence)} base{'' if len(sequence)==1 else 's'}</span>
    #                 </div>
    #                 <div class="no-motif-msg">No motifs detected in this region.</div>
    #                 <div class="sequence-content">
    #                     <span style="color:#555;">{sequence}</span>
    #                 </div>
    #             </div>
    #         """, unsafe_allow_html=True)
    #         return

    #     # Otherwise, proceed as before
    #     ranges = self.parse_motif_range(spans)
    #     highlighted_sequence = ""
    #     previous_end = 0

    #     # Build highlighted sequence
    #     for idx, (start, end) in enumerate(ranges):
    #         motif = motif_ids[idx]
    #         color = motif_colors[int(motif)]
    #         motif_name = motif_names[int(motif)]

    #         # Add interruption if any, use red (#FF0000)
    #         if start > previous_end:
    #             interruption_sequence = sequence[previous_end:start]
    #             highlighted_sequence += (
    #                 f"<span class='interruption' style='color:#fff; background-color:#FF0000; opacity:0.8; font-weight:500;' title='Interruption region ({len(interruption_sequence)} bases)'>"
    #                 f"{interruption_sequence}</span>"
    #             )

    #         # Add motif
    #         motif_sequence = sequence[start:end+1]
    #         highlighted_sequence += (
    #             f"<span class='motif-highlight' style='background-color:{color}' "
    #             f"title='{motif_name} ({len(motif_sequence)} bases)'>"
    #             f"{motif_sequence}</span>"
    #         )
    #         previous_end = end + 1

    #     # Add remaining sequence as interruption, use red (#FF0000)
    #     if previous_end < len(sequence):
    #         interruption_sequence = sequence[previous_end:]
    #         highlighted_sequence += (
    #             f"<span class='interruption' style='color:#fff; background-color:#FF0000; opacity:0.8; font-weight:500;' title='Interruption region ({len(interruption_sequence)} bases)'>"
    #             f"{interruption_sequence}</span>"
    #         )

    #     if sequence_name == "Ref":
    #         sequence_name += "seq"

    #     # visulize the motif sizes 
    #     motif_sizes = [len(motif_names[int(motif)]) for motif in motif_ids]
    #     motif_sizes = sorted(motif_sizes)
    #     motif_sizes = [str(size) for size in motif_sizes]
    #     motif_sizes = ", ".join(motif_sizes)
    #     st.markdown(f"""
    #         <div style="font-size: 12px; color: #555;">
    #             <strong>Motif Sizes:</strong> {motif_sizes}
    #         </div>
    #     """, unsafe_allow_html=True)

    #     st.markdown(f"""
    #         <style>
    #             .sequence-container {{
    #                 font-family: 'SF Mono', 'Monaco', 'Consolas', 'Roboto Mono', monospace;
    #                 margin-bottom: 20px;
    #                 border: 1px solid #e1e5e9;
    #                 border-radius: 12px;
    #                 overflow: hidden;
    #                 box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    #                 background: white;
    #             }}
    #             .sequence-header {{
    #                 background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    #                 color: white;
    #                 padding: 12px 16px;
    #                 font-weight: 600;
    #                 font-size: 14px;
    #                 display: flex;
    #                 justify-content: space-between;
    #                 align-items: center;
    #             }}
    #             .sequence-length {{
    #                 font-size: 12px;
    #                 opacity: 0.9;
    #                 font-weight: normal;
    #             }}
    #             .sequence-content {{
    #                 padding: 16px;
    #                 max-width: 100%;
    #                 white-space: nowrap;
    #                 overflow-x: auto;
    #                 overflow-y: hidden;
    #                 background: #f8f9fa;
    #                 border-top: 1px solid #e1e5e9;
    #                 line-height: 1.6;
    #                 scrollbar-width: thin;  /* Firefox */
    #                 scrollbar-color: #b5b5b5 #f8f9fa; /* Firefox */
    #             }}
    #             .sequence-content::-webkit-scrollbar {{
    #                 height: 10px;
    #             }}
    #             .sequence-content::-webkit-scrollbar-thumb {{
    #                 background: #b5b5b5;
    #                 border-radius: 8px;
    #             }}
    #             .sequence-content::-webkit-scrollbar-track {{
    #                 background: #f8f9fa;
    #             }}
    #             .motif-highlight {{
    #                 display: inline-block;
    #                 padding: 4px 6px;
    #                 margin: 1px;
    #                 border-radius: 6px;
    #                 font-weight: 500;
    #                 color: #1a1a1a;
    #                 box-shadow: 0 1px 3px rgba(0,0,0,0.2);
    #                 transition: all 0.2s ease;
    #                 cursor: pointer;
    #                 border: 1px solid rgba(0,0,0,0.1);
    #             }}
    #             .motif-highlight:hover {{
    #                 transform: translateY(-1px);
    #                 box-shadow: 0 3px 8px rgba(0,0,0,0.3);
    #                 z-index: 10;
    #                 position: relative;
    #             }}
    #             .interruption {{
    #                 display: inline-block;
    #                 padding: 4px 2px;
    #                 margin: 1px;
    #                 color: #fff !important;
    #                 background-color: #FF0000 !important;
    #                 opacity: 0.8;
    #                 font-weight: 500;
    #                 border-radius: 4px;
    #             }}
    #         </style>

    #         <div class="sequence-container">
    #             <div class="sequence-header">
    #                 <span>{sequence_name}</span>
    #                 <span class="sequence-length">Length: {len(sequence)} bases</span>
    #             </div>
    #             <div class="sequence-content">
    #                 {highlighted_sequence}
    #             </div>
    #         </div>
    #     """, unsafe_allow_html=True)


    def display_dynamic_sequence_with_highlighted_motifs(self, sequence_name, sequence, motif_ids, spans, motif_colors, motif_names):
        # Handle the situations where motifs are not available or sequence is just one base
        if (motif_ids == ["."]) or (isinstance(sequence, str) and len(sequence) <= 1):
            if sequence_name == "Ref":
                sequence_name += "seq"
            st.markdown(f"""
                <style>
                    .sequence-dashboard {{
                        font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                        margin-bottom: 25px;
                    }}
                    .sequence-container {{
                        border: 1px solid #e1e5e9;
                        border-radius: 16px;
                        overflow: hidden;
                        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
                        background: white;
                        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                    }}
                    .sequence-container:hover {{
                        box-shadow: 0 8px 30px rgba(0,0,0,0.12);
                        transform: translateY(-2px);
                    }}
                    .sequence-header {{
                        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                        color: white;
                        padding: 16px 20px;
                        font-weight: 700;
                        font-size: 16px;
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                    }}
                    .sequence-length {{
                        font-size: 13px;
                        opacity: 0.95;
                        font-weight: 500;
                        background: rgba(255,255,255,0.15);
                        padding: 4px 12px;
                        border-radius: 20px;
                    }}
                    .no-motif-section {{
                        padding: 20px;
                        background: linear-gradient(135deg, #fff7f7 0%, #feebeb 100%);
                        border: 1px solid #feb2b2;
                        border-radius: 12px;
                        margin: 20px;
                        text-align: center;
                    }}
                    .no-motif-icon {{
                        font-size: 24px;
                        margin-bottom: 8px;
                    }}
                    .no-motif-msg {{
                        color: #c53030;
                        font-weight: 600;
                        font-size: 15px;
                        margin-bottom: 8px;
                    }}
                    .no-motif-desc {{
                        color: #744210;
                        font-size: 13px;
                        opacity: 0.8;
                    }}
                    .sequence-content {{
                        padding: 20px;
                        max-width: 100%;
                        white-space: nowrap;
                        overflow-x: auto;
                        overflow-y: hidden;
                        background: #f8fafc;
                        line-height: 1.8;
                        scrollbar-width: thin;
                        scrollbar-color: #cbd5e0 #f8fafc;
                        font-size: 15px;
                        font-weight: 500;
                    }}
                    .sequence-content::-webkit-scrollbar {{
                        height: 8px;
                    }}
                    .sequence-content::-webkit-scrollbar-thumb {{
                        background: #cbd5e0;
                        border-radius: 10px;
                    }}
                    .sequence-content::-webkit-scrollbar-track {{
                        background: #f8fafc;
                    }}
                </style>
                
                <div class="sequence-dashboard">
                    <div class="sequence-container">
                        <div class="sequence-header">
                            <span>{sequence_name}</span>
                            <span class="sequence-length">{len(sequence)} base{'' if len(sequence)==1 else 's'}</span>
                        </div>
                        <div class="no-motif-section">
                            <div class="no-motif-icon">🔍</div>
                            <div class="no-motif-msg">No motifs detected in this region</div>
                            <div class="no-motif-desc">The sequence contains no recognizable motif patterns</div>
                        </div>
                        <div class="sequence-content">
                            <span style="color:#4a5568; letter-spacing:0.5px;">{sequence}</span>
                        </div>
                    </div>
                </div>
            """, unsafe_allow_html=True)
            return

        # Process sequences with motifs
        ranges = self.parse_motif_range(spans)
        highlighted_sequence = ""
        previous_end = 0

        # Create legend
        legend_html = ""
        unique_motifs = sorted(set(motif_ids))
        for motif_id in unique_motifs:
            color = motif_colors[int(motif_id)]
            motif_name = motif_names[int(motif_id)]
            legend_html += f"""
            <div class="legend-item" data-motif="{motif_id}">
                <span class="legend-color" style="background-color:{color};"></span>
                <span class="legend-text">{motif_name}</span>
            </div>
            """

        # Build highlighted sequence with interactive elements
        for idx, (start, end) in enumerate(ranges):
            motif = motif_ids[idx]
            color = motif_colors[int(motif)]
            motif_name = motif_names[int(motif)]

            # Add interruption if any
            if start > previous_end:
                interruption_sequence = sequence[previous_end:start]
                highlighted_sequence += (
                    f"<span class='interruption-segment' data-type='interruption' "
                    f"data-content='Interruption Region | Length: {len(interruption_sequence)} bases | Sequence: {interruption_sequence}'>"
                    f"<span class='interruption-text'>{interruption_sequence}</span>"
                    f"</span>"
                )

            # Add motif
            motif_sequence = sequence[start:end+1]
            highlighted_sequence += (
                f"<span class='motif-segment motif-{motif}' data-motif='{motif}' "
                f"style='background-color:{color};' "
                f"data-content='{motif_name} | Length: {len(motif_sequence)} bases | Sequence: {motif_sequence}'>"
                f"<span class='motif-text'>{motif_sequence}</span>"
                f"<span class='motif-badge'>{motif_name}</span>"
                f"</span>"
            )
            previous_end = end + 1

        # Add remaining sequence as interruption
        if previous_end < len(sequence):
            interruption_sequence = sequence[previous_end:]
            highlighted_sequence += (
                f"<span class='interruption-segment' data-type='interruption' "
                f"data-content='Interruption Region | Length: {len(interruption_sequence)} bases | Sequence: {interruption_sequence}'>"
                f"<span class='interruption-text'>{interruption_sequence}</span>"
                f"</span>"
            )

        if sequence_name == "Ref":
            sequence_name += "seq"

        # Calculate motif statistics
        # Collect motif sizes in the order motifs appear in the legend (i.e., motif_names order, which matches the legend)
        legend_motif_sizes = [len(name) for name in motif_names]
        # Remove duplicates while preserving order
        seen_sizes = set()
        ordered_unique_sizes = []
        for size in legend_motif_sizes:
            if size not in seen_sizes:
                ordered_unique_sizes.append(size)
                seen_sizes.add(size)
        motif_sizes_display = ", ".join([f"{size} bp" for size in ordered_unique_sizes])
        total_motifs = len(motif_ids)
        coverage = sum([end - start + 1 for start, end in ranges]) / len(sequence) * 100

        st.html(f"""
            <style>
                .sequence-dashboard {{
                    font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                    margin-bottom: 25px;
                }}
                .sequence-container {{
                    border: 1px solid #e1e5e9;
                    border-radius: 16px;
                    overflow: hidden;
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08);
                    background: white;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                }}
                .sequence-container:hover {{
                    box-shadow: 0 8px 30px rgba(0,0,0,0.12);
                }}
                .sequence-header {{
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                    padding: 16px 20px;
                    font-weight: 700;
                    font-size: 16px;
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                }}
                .sequence-length {{
                    font-size: 13px;
                    opacity: 0.95;
                    font-weight: 500;
                    background: rgba(255,255,255,0.15);
                    padding: 4px 12px;
                    border-radius: 20px;
                }}
                .stats-bar {{
                    display: flex;
                    gap: 20px;
                    padding: 15px 20px;
                    background: #f0f4ff;
                    border-bottom: 1px solid #e1e8ff;
                    font-size: 13px;
                    color: #4a5568;
                }}
                .stat-item {{
                    display: flex;
                    align-items: center;
                    gap: 6px;
                    font-weight: 500;
                }}
                .stat-value {{
                    background: white;
                    padding: 2px 8px;
                    border-radius: 12px;
                    border: 1px solid #cbd5e0;
                    font-weight: 600;
                    color: #2d3748;
                }}
                .legend-container {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 12px;
                    padding: 15px 20px;
                    background: #f8fafc;
                    border-bottom: 1px solid #e2e8f0;
                }}
                .legend-title {{
                    font-weight: 600;
                    color: #4a5568;
                    font-size: 13px;
                    display: flex;
                    align-items: center;
                }}
                .legend-item {{
                    display: flex;
                    align-items: center;
                    gap: 8px;
                    padding: 6px 12px;
                    border-radius: 8px;
                    background: white;
                    border: 2px solid transparent;
                    cursor: pointer;
                    transition: all 0.2s ease;
                    font-size: 13px;
                    font-weight: 500;
                }}
                .legend-item:hover {{
                    transform: translateY(-1px);
                    box-shadow: 0 4px 12px rgba(0,0,0,0.1);
                    border-color: #667eea;
                }}
                .legend-color {{
                    width: 16px;
                    height: 16px;
                    border-radius: 4px;
                    border: 2px solid white;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                .sequence-content {{
                    padding: 20px;
                    max-width: 100%;
                    white-space: nowrap;
                    overflow-x: auto;
                    overflow-y: hidden;
                    background: #f8fafc;
                    line-height: 1.8;
                    scrollbar-width: thin;
                    scrollbar-color: #cbd5e0 #f8fafc;
                    font-size: 15px;
                    font-weight: 500;
                    min-height: 80px;
                }}
                .sequence-content::-webkit-scrollbar {{
                    height: 8px;
                }}
                .sequence-content::-webkit-scrollbar-thumb {{
                    background: #cbd5e0;
                    border-radius: 10px;
                }}
                .sequence-content::-webkit-scrollbar-track {{
                    background: #f8fafc;
                }}
                .motif-segment {{
                    display: inline-flex;
                    align-items: center;
                    gap: 4px;
                    padding: 6px 8px;
                    margin: 2px;
                    border-radius: 10px;
                    font-weight: 600;
                    color: #1a202c;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.15);
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                    cursor: pointer;
                    border: 2px solid rgba(255,255,255,0.3);
                    position: relative;
                    overflow: hidden;
                }}
                .motif-segment:hover {{
                    transform: translateY(-2px) scale(1.02);
                    box-shadow: 0 8px 25px rgba(0,0,0,0.25);
                    z-index: 10;
                }}
                .motif-segment::before {{
                    content: '';
                    position: absolute;
                    top: 0;
                    left: -100%;
                    width: 100%;
                    height: 100%;
                    background: linear-gradient(90deg, transparent, rgba(255,255,255,0.4), transparent);
                    transition: left 0.6s ease;
                }}
                .motif-segment:hover::before {{
                    left: 100%;
                }}
                .motif-text {{
                    letter-spacing: 0.5px;
                }}
                .interruption-segment {{
                    display: inline-block;
                    padding: 6px 4px;
                    margin: 2px;
                    color: #fff;
                    background: linear-gradient(135deg, #fc8181, #e53e3e);
                    opacity: 0.9;
                    font-weight: 500;
                    border-radius: 6px;
                    cursor: pointer;
                    transition: all 0.3s ease;
                    border: 1px dashed rgba(255,255,255,0.5);
                }}
                .interruption-segment:hover {{
                    opacity: 1;
                    transform: scale(1.05);
                    box-shadow: 0 4px 12px rgba(229, 62, 62, 0.3);
                }}
                .interruption-text {{
                    letter-spacing: 0.5px;
                }}
                .tooltip {{
                    position: fixed;
                    background: rgba(45, 55, 72, 0.95);
                    color: white;
                    padding: 10px 14px;
                    border-radius: 10px;
                    font-size: 12px;
                    pointer-events: none;
                    opacity: 0;
                    transform: translateY(10px);
                    transition: all 0.2s ease;
                    z-index: 1000;
                    backdrop-filter: blur(10px);
                    border: 1px solid rgba(255,255,255,0.1);
                    max-width: 300px;
                    font-family: 'SF Mono', monospace;
                }}
                .tooltip.show {{
                    opacity: 1;
                    transform: translateY(0);
                }}
                .highlighted {{
                    transform: scale(1.1);
                    box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.5);
                    z-index: 5;
                }}
            </style>

            <div class="sequence-dashboard">
                <div class="sequence-container">
                    <div class="sequence-header">
                        <span>{sequence_name}</span>
                        <span class="sequence-length">{len(sequence)} bases</span>
                    </div>
                    
                    <div class="stats-bar">
                        <div class="stat-item">
                            <span>Motifs:</span>
                            <span class="stat-value">{total_motifs}</span>
                        </div>
                        <div class="stat-item">
                            <span>Coverage:</span>
                            <span class="stat-value">{coverage:.1f}%</span>
                        </div>
                        <div class="stat-item">
                            <span>Motif Sizes:</span>
                            <span class="stat-value">{motif_sizes_display}</span>
                        </div>
                    </div>
                    
                    <div class="legend-container">
                        <div class="legend-title">Motifs Legend:</div>
                        {legend_html}
                        <div class="legend-item">
                            <div class="legend-color" style="background:linear-gradient(135deg, #fc8181, #e53e3e);"></div>
                            <span class="legend-text">Interruption</span>
                        </div>
                    </div>

                    <div class="sequence-content" id="sequence-content-{sequence_name}">
                        {highlighted_sequence}
                    </div>
                </div>
            </div>

            <script>
                function initializeSequenceVisualization(seqName) {{
                    const content = document.getElementById('sequence-content-' + seqName);
                    const tooltip = document.createElement('div');
                    tooltip.className = 'tooltip';
                    document.body.appendChild(tooltip);
                    
                    const legendItems = document.querySelectorAll('.legend-item[data-motif]');
                    
                    function showTooltip(event) {{
                        const target = event.target.closest('[data-content]');
                        if (target && target.dataset.content) {{
                            tooltip.innerHTML = target.dataset.content.replace(/\|/g, '<br>');
                            tooltip.style.left = (event.pageX + 15) + 'px';
                            tooltip.style.top = (event.pageY + 15) + 'px';
                            tooltip.classList.add('show');
                        }}
                    }}
                    
                    function hideTooltip() {{
                        tooltip.classList.remove('show');
                    }}
                    
                    function highlightMotif(motifId) {{
                        document.querySelectorAll('.motif-segment').forEach(segment => {{
                            if (segment.classList.contains('motif-' + motifId)) {{
                                segment.classList.add('highlighted');
                            }} else {{
                                segment.style.opacity = '0.5';
                            }}
                        }});
                    }}
                    
                    function resetHighlight() {{
                        document.querySelectorAll('.motif-segment').forEach(segment => {{
                            segment.classList.remove('highlighted');
                            segment.style.opacity = '1';
                        }});
                    }}
                    
                    content.addEventListener('mousemove', showTooltip);
                    content.addEventListener('mouseleave', hideTooltip);
                    
                    legendItems.forEach(item => {{
                        item.addEventListener('mouseover', () => highlightMotif(item.dataset.motif));
                        item.addEventListener('mouseout', resetHighlight);
                    }});
                }}
                
                initializeSequenceVisualization('{sequence_name}');
            </script>
        """)
    
    # def display_dynamic_sequence_with_highlighted_motifs(self, sequence_name, sequence, motif_ids, spans, motif_colors, motif_names):

    #     ranges = self.parse_motif_range(spans)
    #     highlighted_sequence = ""
    #     previous_end = 0
    #     for idx, (start, end) in enumerate(ranges):
    #         motif = motif_ids[idx]
    #         color = motif_colors[int(motif)]
    #         motif_name = motif_names[int(motif)]
    #         if start > previous_end:
    #             interruption_sequence = sequence[previous_end:start]
    #             highlighted_sequence += (
    #                 f"<span style='background-color:#FF0000; padding:2px; border-radius:4px;' title='Interruption'>"
    #                 f"{interruption_sequence}</span>"
    #             )
            
    #         motif_sequence = sequence[start:end+1]
    #         highlighted_sequence += (
    #             f"<span style='background-color:{color}; padding:2px; border-radius:4px;' title='Motif: {motif_name}'>"
    #             f"{motif_sequence}</span>"
    #         )

    #         previous_end = end + 1

    #     if previous_end < len(sequence):
    #         interruption_sequence = sequence[previous_end:]
    #         highlighted_sequence += (
    #             f"<span style='background-color:#FF0000; padding:2px; border-radius:4px;'title='Interruption'>"
    #             f"{interruption_sequence}</span>"
    #         )
    #     if sequence_name == "Ref":
    #         sequence_name += "seq"
    #     st.markdown(f"""
    #         <div style="display: flex; align-items: right;">
    #             <div style="font-family:monospace; font-size:16px; padding:10px; border:1px solid black; border-radius:8px; margin-right: 10px; text-align: right;">
    #                 <strong>{sequence_name}</strong>
    #             </div>
    #             <div style="font-family:monospace; font-size:16px; width:100%; max-height:120px; overflow-x:auto; white-space:nowrap; padding:10px; border:1px solid black; border-radius:8px;">
    #                 {highlighted_sequence}
    #             </div>
    #         </div>
    #     """, unsafe_allow_html=True)


    
    # def display_motifs_as_bars(self,sequence_name, motif_colors, motif_ids, spans, sequence, motif_names):
    #     sequence_length = len(sequence)
    #     ranges = self.parse_motif_range(spans)
        
    #     if not isinstance(motif_names, list):
    #         motif_names = [motif_names]
        
    #     bar_container = "<div style='width:100%; position: relative; height: 30px; border:2px solid black; border-radius: 8px;'>"
    #     previous_end = 0
    #     gap = 0.05  # Small gap between motifs for better visibility

    #     for idx, (start, end) in enumerate(ranges):
    #         motif = motif_ids[idx]
    #         color = motif_colors[int(motif)]
    #         span_length = end - start + 1  

    #         if start >= 0 and end <= sequence_length:
    #             if start > previous_end:
    #                 interruption_width = (start - previous_end) / sequence_length * 100
    #                 interruption_start = previous_end / sequence_length * 100
    #                 bar_container += (
    #                     f"<div style='position:absolute; background-color:#FF0000; left:{interruption_start}%; "
    #                     f"width:{interruption_width}%; height:28px; top:-1px; border-radius:6px; border:1px solid black; "
    #                     f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2); cursor: pointer;' "
    #                     f"title='Interruption: {sequence[previous_end:start]}'>"
    #                     f"</div>"
    #                 )

    #             relative_width = (span_length / sequence_length) * 100 - gap
    #             relative_start = (start / sequence_length) * 100

    #             bar_container += (
    #                 f"<div class='hoverable-div' style='position:absolute; background-color:{color}; left:{relative_start}%; "
    #                 f"width:{relative_width}%; height:28px; top:-1px; border-radius:6px; border:1px solid black; "
    #                 f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2); cursor: pointer;' "
    #                 f"title='Motif: {motif_names[int(motif)]}'>"
    #                 f"</div>"
    #             )

    #             previous_end = end + 1

    #     # Add interruption bar after the last motif
    #     if previous_end < sequence_length:
    #         interruption_width = (sequence_length - previous_end) / sequence_length * 100
    #         interruption_start = previous_end / sequence_length * 100
    #         bar_container += (
    #             f"<div style='position:absolute; background-color:#FF0000; left:{interruption_start}%; "
    #             f"width:{interruption_width}%; height:28px; top:-1px; border-radius:6px; border:1px solid black; "
    #             f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2); cursor: pointer;' "
    #             f"title='Interruption: {sequence[previous_end:]}'>"
    #             f"</div>"
    #         )

    #     bar_container += "</div>"
    #     if sequence_name == "Ref":
    #         sequence_name += "seq"

    #     st.markdown(f"""
    #         <div style="display: flex; align-items: center;">
    #             <div style="font-family:monospace; font-size:16px; padding:10px; border:1px solid black; border-radius:8px; margin-right: 10px;">
    #                 <strong>{sequence_name}</strong>
    #             </div>
    #             {bar_container}
    #         </div>
    #     """, unsafe_allow_html=True)

    def display_motifs_as_bars(self, sequence_name, motif_colors, motif_ids, spans, sequence, motif_names):
        sequence_length = len(sequence)
        ranges = self.parse_motif_range(spans)
        
        if not isinstance(motif_names, list):
            motif_names = [motif_names]
        
        # Create legend
        legend_html = ""
        unique_motifs = sorted(set(motif_ids))
        for motif_id in unique_motifs:
            color = motif_colors[int(motif_id)]
            motif_name = motif_names[int(motif_id)]
            legend_html += f"""
            <div class="legend-item" data-motif="{motif_id}">
                <span class="legend-color" style="background-color:{color};"></span>
                <span class="legend-text">{motif_name}</span>
            </div>
            """

        bar_container = f"""
        <div class="bar-visualization">
            <div class="sequence-bar-container" id="bar-container-{sequence_name}">
                <div class="bar-track">
        """
        
        previous_end = 0
        gap = 0.3

        for idx, (start, end) in enumerate(ranges):
            motif = motif_ids[idx]
            color = motif_colors[int(motif)]
            span_length = end - start + 1  
            motif_name = motif_names[int(motif)]

            if start >= 0 and end <= sequence_length:
                if start > previous_end:
                    interruption_width = (start - previous_end) / sequence_length * 100
                    interruption_start = previous_end / sequence_length * 100
                    bar_container += (
                        f"<div class='interruption-bar' data-motif='interruption' "
                        f"style='left:{interruption_start}%; width:{interruption_width}%;' "
                        f"data-content='Interruption: {sequence[previous_end:start]} ({start-previous_end} bases)'>"
                        f"<div class='bar-pattern'></div>"
                        f"</div>"
                    )

                relative_width = max((span_length / sequence_length) * 100 - gap, 0.5)
                relative_start = (start / sequence_length) * 100

                bar_container += (
                    f"<div class='motif-bar motif-{motif}' data-motif='{motif}' "
                    f"style='left:{relative_start}%; width:{relative_width}%; background-color:{color};' "
                    f"data-content='{motif_name}: {sequence[start:end+1]} ({span_length} bases)'>"
                    f"<div class='bar-glow'></div>"
                    f"<div class='bar-label'>{motif_name if motif_name else 'M'}</div>"
                    f"</div>"
                )

                previous_end = end + 1

        if previous_end < sequence_length:
            interruption_width = (sequence_length - previous_end) / sequence_length * 100
            interruption_start = previous_end / sequence_length * 100
            bar_container += (
                f"<div class='interruption-bar' data-motif='interruption' "
                f"style='left:{interruption_start}%; width:{interruption_width}%;' "
                f"data-content='Interruption: {sequence[previous_end:]} ({sequence_length-previous_end} bases)'>"
                f"<div class='bar-pattern'></div>"
                f"</div>"
            )

        bar_container += """
                </div>
                <div class="bar-scale">
                    <span>0</span>
                    <span style="left:25%">25%</span>
                    <span style="left:50%">50%</span>
                    <span style="left:75%">75%</span>
                    <span style="right:0">100%</span>
                </div>
            </div>
            <div id="tooltip-""" + sequence_name + """\" class="tooltip"></div>
        </div>
        """

        if sequence_name == "Ref":
            sequence_name += "seq"

        st.html(f"""
            <style>
                .bar-visualization {{
                    font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                    margin: 20px 0;
                }}
                .sequence-header {{
                    display: flex;
                    align-items: center;
                    justify-content: space-between;
                    margin-bottom: 15px;
                    padding: 0 10px;
                }}
                .sequence-title {{
                    font-size: 18px;
                    font-weight: 700;
                    color: #2d3748;
                    background: linear-gradient(135deg, #667eea, #764ba2);
                    -webkit-background-clip: text;
                    -webkit-text-fill-color: transparent;
                    padding: 8px 16px;
                    border-radius: 12px;
                    border: 2px solid #e2e8f0;
                }}
                .sequence-info {{
                    font-size: 14px;
                    color: #718096;
                    font-weight: 500;
                }}
                .sequence-bar-container {{
                    position: relative;
                    height: 60px;
                    background: #f7fafc;
                    border-radius: 16px;
                    border: 2px solid #e2e8f0;
                    padding: 10px;
                    margin: 15px 0;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                }}
                .sequence-bar-container:hover {{
                    border-color: #667eea;
                    box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15);
                    transform: translateY(-2px);
                }}
                .bar-track {{
                    position: relative;
                    height: 40px;
                    background: #edf2f7;
                    border-radius: 12px;
                    overflow: hidden;
                    border: 1px solid #cbd5e0;
                }}
                .motif-bar {{
                    position: absolute;
                    height: 36px;
                    top: 2px;
                    border-radius: 10px;
                    border: 2px solid white;
                    box-shadow: 0 4px 12px rgba(0,0,0,0.1);
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                    cursor: pointer;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    font-weight: 700;
                    font-size: 12px;
                    color: white;
                    text-shadow: 0 1px 2px rgba(0,0,0,0.3);
                    overflow: hidden;
                }}
                .motif-bar:hover {{
                    transform: scale(1.05) translateY(-2px);
                    box-shadow: 0 8px 25px rgba(0,0,0,0.25);
                    z-index: 10;
                }}
                .bar-glow {{
                    position: absolute;
                    top: 0;
                    left: -100%;
                    width: 100%;
                    height: 100%;
                    background: linear-gradient(90deg, transparent, rgba(255,255,255,0.4), transparent);
                    transition: left 0.6s ease;
                }}
                .motif-bar:hover .bar-glow {{
                    left: 100%;
                }}
                .interruption-bar {{
                    position: absolute;
                    height: 36px;
                    top: 2px;
                    background: #fed7d7;
                    border: 2px dashed #e53e3e;
                    border-radius: 8px;
                    cursor: pointer;
                    transition: all 0.3s ease;
                    overflow: hidden;
                }}
                .interruption-bar:hover {{
                    background: #feebeb;
                    transform: scaleY(1.1);
                }}
                .bar-pattern {{
                    width: 100%;
                    height: 100%;
                    background: repeating-linear-gradient(
                        45deg,
                        transparent,
                        transparent 5px,
                        rgba(229, 62, 62, 0.1) 5px,
                        rgba(229, 62, 62, 0.1) 10px
                    );
                }}
                .bar-scale {{
                    position: relative;
                    height: 20px;
                    margin-top: 10px;
                    font-size: 11px;
                    color: #718096;
                    font-weight: 500;
                }}
                .bar-scale span {{
                    position: absolute;
                    transform: translateX(-50%);
                }}
                .bar-scale span:first-child {{ left: 0; transform: none; }}
                .bar-scale span:last-child {{ left: auto; right: 0; transform: none; }}
                .tooltip {{
                    position: fixed;
                    background: rgba(45, 55, 72, 0.95);
                    color: white;
                    padding: 8px 12px;
                    border-radius: 8px;
                    font-size: 12px;
                    pointer-events: none;
                    opacity: 0;
                    transform: translateY(10px);
                    transition: all 0.2s ease;
                    z-index: 1000;
                    backdrop-filter: blur(10px);
                    border: 1px solid rgba(255,255,255,0.1);
                    max-width: 300px;
                }}
                .tooltip.show {{
                    opacity: 1;
                    transform: translateY(0);
                }}
                .legend-container {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 12px;
                    margin: 15px 0;
                    padding: 15px;
                    background: #f7fafc;
                    border-radius: 12px;
                    border: 1px solid #e2e8f0;
                }}
                .legend-item {{
                    display: flex;
                    align-items: center;
                    gap: 8px;
                    padding: 6px 12px;
                    border-radius: 8px;
                    background: white;
                    border: 2px solid transparent;
                    cursor: pointer;
                    transition: all 0.2s ease;
                    font-size: 13px;
                    font-weight: 500;
                }}
                .legend-item:hover {{
                    transform: translateY(-2px);
                    box-shadow: 0 4px 12px rgba(0,0,0,0.1);
                    border-color: #667eea;
                }}
                .legend-color {{
                    width: 16px;
                    height: 16px;
                    border-radius: 4px;
                    border: 2px solid white;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                .highlighted {{
                    transform: scale(1.1);
                    box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.3);
                    z-index: 5;
                }}
            </style>

            <div class="bar-visualization">
                <div class="sequence-header">
                    <div class="sequence-title">{sequence_name}</div>
                    <div class="sequence-info">{sequence_length} bases • {len(ranges)} motifs</div>
                </div>
                
                <div class="legend-container">
                    <div style="font-weight:600; color:#4a5568;">Legend:</div>
                    {legend_html}
                    <div class="legend-item">
                        <div class="legend-color" style="background:#fed7d7; border:2px dashed #e53e3e;"></div>
                        <span class="legend-text">Interruption</span>
                    </div>
                </div>

                {bar_container}
            </div>

            <script>
                function initializeBarVisualization(sequenceName) {{
                    const container = document.getElementById('bar-container-' + sequenceName);
                    const tooltip = document.getElementById('tooltip-' + sequenceName);
                    const legendItems = document.querySelectorAll('.legend-item[data-motif]');
                    
                    function showTooltip(event) {{
                        const target = event.target.closest('[data-content]');
                        if (target && target.dataset.content) {{
                            tooltip.innerHTML = target.dataset.content;
                            tooltip.style.left = (event.pageX + 10) + 'px';
                            tooltip.style.top = (event.pageY + 10) + 'px';
                            tooltip.classList.add('show');
                        }}
                    }}
                    
                    function hideTooltip() {{
                        tooltip.classList.remove('show');
                    }}
                    
                    function highlightMotif(motifId) {{
                        document.querySelectorAll('.motif-bar').forEach(bar => {{
                            if (bar.classList.contains('motif-' + motifId)) {{
                                bar.classList.add('highlighted');
                            }} else {{
                                bar.style.opacity = '0.4';
                            }}
                        }});
                    }}
                    
                    function resetHighlight() {{
                        document.querySelectorAll('.motif-bar').forEach(bar => {{
                            bar.classList.remove('highlighted');
                            bar.style.opacity = '1';
                        }});
                    }}
                    
                    container.addEventListener('mousemove', showTooltip);
                    container.addEventListener('mouseleave', hideTooltip);
                    
                    legendItems.forEach(item => {{
                        item.addEventListener('mouseover', () => highlightMotif(item.dataset.motif));
                        item.addEventListener('mouseout', resetHighlight);
                    }});
                }}
                
                initializeBarVisualization('{sequence_name}');
            </script>
        """)
        
    def display_motifs_with_bars(self, record, left_column, right_column, motif_colors, CN1_col, CN2_col, show_comparison, hgsvc_records):
        motif_names, motif_count_h1, motif_count_h2 = self.parse_motif_in_region(record)
        if motif_count_h1 == None and motif_count_h2 == None:
            st.info("No motifs found in the region")
            return
        CN1_col.markdown(f""" 
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 1 Total copy number:</strong> {str(record['spans'][1]).count('-')}
            </div>
        """, unsafe_allow_html=True)

        with right_column:
            
            self.display_motif_legend(motif_names, motif_colors, right_column)

        if record['alt_allele2'] != '':
            CN2_col.markdown(f"""
                <div style="font-size: 20px; color: #FF5733;">
                    <strong>Allele 2 Total copy number:</strong> {str(record['spans'][2]).count('-')}
                </div>
            """, unsafe_allow_html=True)


        with left_column:
            tab1, tab2 , tab3 = st.tabs(["Alleles", "Alleles vs Ref", "Alleles vs Pop"])
            st.markdown(
                "<style>.tab-content {font-size: 20px;}</style>",
                unsafe_allow_html=True,
            )
            with tab2:
                self.display_motifs_as_bars("Ref", motif_colors, record['motif_ids_ref'], record['spans'][0], record['ref_allele'], motif_names)
                self.display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names)
                if record['alt_allele2'] != '':
                    self.display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names)
            with tab1:
                self.display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names)
                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        plot_motif_bar(motif_count_h1, motif_names, motif_colors)
                
                if record['alt_allele2'] != '':
                    self.display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names,)
                    plot_container_h2 = st.empty()

                    with plot_container_h2:
                        if show_comparison == False:
                            plot_motif_bar(motif_count_h2, motif_names, motif_colors)

            with tab3:
                if hgsvc_records:
                    self.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
                else:
                    st.info("no population data found")

            

    # def plot_motif_bar(self, motif_count, motif_names, motif_colors=None):
    #     motif_labels = []
    #     motif_counts = []
    #     for label, value in sorted(motif_count.items()):
    #         motif_name = motif_names[int(label)]
    #         if motif_name:
    #             motif_labels.append(motif_name)
    #             motif_counts.append(value)
            
    #     data = {
    #         'Motif': motif_labels,
    #         'Count': motif_counts
    #     }
    #     df = pd.DataFrame(data)
        
    #     color_list = [motif_colors[int(label)] for label in sorted(motif_count.keys()) if motif_colors is not None and int(label) < len(motif_colors)]    
    #     bar_chart = alt.Chart(df).mark_bar().encode(
    #         x=alt.X('Motif', sort=None),
    #         y='Count',
    #         tooltip=['Motif', 'Count'],
    #         color=alt.Color('Motif', scale=alt.Scale(domain=motif_labels, range=color_list))  # Ensure the scale matches correctly
    #     ).properties(
    #         width=200,
    #         height=200,
    #         title="Motif Occurrences"
    #     )

    #     bar_chart = bar_chart.configure_axisX(labelAngle=0)
    #     st.altair_chart(bar_chart, use_container_width=True)

    def plot_Cohort_results(self,cohort_records):
        sequences = []
        span_list = []
        motif_ids_list = []
        # make space between the last print 
        sort_by = st.radio("Sort by:", ("Value", "Sample Name"), horizontal=True, key="sort_by_cohort")
        for key in cohort_records.keys():
            sequences.append({'name': f'{key}_alle1', 'sequence': cohort_records[key]['alt_allele1']})
            span_list.append(cohort_records[key]['spans'][1])
            motif_ids_list.append(cohort_records[key]['motif_ids_h1'])
            if cohort_records[key]['alt_allele2'] != '':
                sequences.append({'name': f'{key}_alle2', 'sequence': cohort_records[key]['alt_allele2']})
                span_list.append(cohort_records[key]['spans'][2])
                motif_ids_list.append(cohort_records[key]['motif_ids_h2'])
            else :
                sequences.append({'name': f'{key}_alle2', 'sequence': ""})
                span_list.append("")
                motif_ids_list.append([0])

        motif_names = cohort_records[list(cohort_records.keys())[0]]['motifs']
        record = cohort_records[list(cohort_records.keys())[0]]
        motif_colors, df = self.stack_plot(record, motif_names, sequences, span_list, motif_ids_list,sort_by)

        
        region = f"{record['chr']}:{record['pos']-1}-{record['stop']-1}"
        self.bar_plot_motif_count(df, region, sort_by=sort_by)

        motif_counts = df[df['Motif'] != 'Interruption'].groupby(['Sample', 'Motif']).size().reset_index(name='Count')

        heatmap_data = motif_counts.pivot(index='Sample', columns='Motif', values='Count').fillna(0)

        heatmap_data_long = heatmap_data.reset_index().melt(id_vars='Sample', var_name='Motif', value_name='Count')

        self.plot_heatmap(heatmap_data_long, sort_by)

        figure = go.Figure()

        if df.empty:
            st.warning("No motifs found in the region")
            st.stop()
        unique_samples = df['Sample'].unique()
        unique_samples = [sample for sample in unique_samples if sample != "Interruption"]
        
        ref_data = []
        allele_data = []


        for sample in unique_samples:
            sample_df = df[df['Sample'] == sample]
            unique_motifs = sample_df['Motif'].unique()
            unique_motifs = [motif for motif in unique_motifs if motif != "Interruption"]
            for motif in unique_motifs:
                motif_df = sample_df[sample_df['Motif'] == motif]
                trace = go.Scatter(
                    x=[motif],
                    y=[len(motif_df)],
                    mode='markers',
                    name=sample,
                    marker=dict(size=20)
                )
                if sample in ["Ref", "Allel1", "Allel2"]:
                    if sample == "Ref":
                        # change the marker to X
                        trace['marker']['symbol'] = "x"
                        ref_data.append(trace)
                    else:
                        trace['marker']['symbol'] = "triangle-down"
                        
                        allele_data.append(trace)
                else:
                    figure.add_trace(trace)

        for trace in ref_data + allele_data:
            figure.add_trace(trace)
        color_mapping = { sample: np.random.choice(px.colors.qualitative.Plotly) for sample in unique_samples}
        color_mapping["Ref"] = "green"

        figure.for_each_trace(lambda trace: trace.update(marker=dict(color=color_mapping.get(trace.name, 'gray'))))

        unique_legend_names = set()
        for trace in figure.data:
            if trace.name in unique_legend_names:
                trace.showlegend = False
            elif trace.name in color_mapping:
                unique_legend_names.add(trace.name)
                trace.showlegend = True
            else:
                trace.showlegend = False
        
        figure.update_layout(
            title="Motif Occurrences",
            xaxis_title="Motif",
            yaxis_title="Count",
        )

        xaxis_colors = {motif: motif_colors[idx] for idx, motif in enumerate(motif_names)}
        figure.update_xaxes(tickmode='array', tickvals=list(xaxis_colors.keys()), ticktext=[
            f'<span style="color:{xaxis_colors[motif]}">{motif}</span>' for motif in xaxis_colors.keys()
        ], tickangle=45)
        figure.update_yaxes(range=[0, df['Sample'].value_counts().max()])

        st.plotly_chart(figure, use_container_width=True)
    def plot_HGSVC_VS_allele(self, record, hgsvc_records, motif_names):
        sequences = []
        span_list = []
        motif_ids_list = []
        sort_by = st.radio("Sort by:", ("Value", "Sample Name"), horizontal=True, key="sort_by_HGSVC")
        sequences.append({'name': "Ref", 'sequence': record['ref_allele']})
        span_list.append(record['spans'][0])
        motif_ids_list.append(record['motif_ids_ref'])
        sequences.append({'name': "Allel1", 'sequence':  record['alt_allele1']})
        span_list.append(record['spans'][1])
        motif_ids_list.append(record['motif_ids_h1'])
        if record['alt_allele2'] != '':
            sequences.append({'name': "Allel2", 'sequence': record['alt_allele2']})
            span_list.append(record['spans'][2])
            motif_ids_list.append(record['motif_ids_h2'])

        for key in hgsvc_records.keys():
            if hgsvc_records[key]['alt_allele'] == '':
                continue
            sequences.append({'name': key, 'sequence': hgsvc_records[key]['alt_allele']})
            span_list.append(hgsvc_records[key]['spans'][1])
            motif_ids_list.append(hgsvc_records[key]['motif_ids_h'])

        motif_colors, df = self.stack_plot(record, motif_names, sequences, span_list, motif_ids_list,sort_by)
        region = f"{record['chr']}:{record['pos']-1}-{record['stop']-1}"
        self.bar_plot_motif_count(df, region, sort_by)




        pivot_hgsvc = pd.pivot_table(df[df['Sample'] == 'HGSVC'], index='Motif', columns='Sample', values='Length', aggfunc='count', fill_value=0)
        pivot_sample = pd.pivot_table(df[df['Sample'] != 'HGSVC'], index='Motif', columns='Sample', values='Length', aggfunc='count', fill_value=0)

        pivot_hgsvc_long = pivot_hgsvc.reset_index().melt(id_vars='Motif', var_name='Sample', value_name='Count')
        pivot_sample_long = pivot_sample.reset_index().melt(id_vars='Motif', var_name='Sample', value_name='Count')

        combined_data = pd.concat([pivot_hgsvc_long, pivot_sample_long])

        self.plot_heatmap(combined_data, sort_by=sort_by)
        
        figure = go.Figure()

        unique_samples = df['Sample'].unique()
        unique_samples = [sample for sample in unique_samples if sample != "Interruption"]

        ref_data = []
        allele_data = []

        for sample in unique_samples:
            sample_df = df[df['Sample'] == sample]
            unique_motifs = sample_df['Motif'].unique()
            unique_motifs = [motif for motif in unique_motifs if motif != "Interruption"]
            for motif in unique_motifs:
                motif_df = sample_df[sample_df['Motif'] == motif]
                trace = go.Scatter(
                    x=[motif],
                    y=[len(motif_df)],
                    mode='markers',
                    name=sample,
                    marker=dict(size=20)
                )
                if sample in ["Ref", "Allel1", "Allel2"]:
                    if sample == "Ref":
                        trace['marker']['symbol'] = "x"
                        ref_data.append(trace)
                    else:
                        trace['marker']['symbol'] = "triangle-down"
                        
                        allele_data.append(trace)
                else:
                    figure.add_trace(trace)

        for trace in ref_data + allele_data:
            figure.add_trace(trace)

        color_mapping = {
            "Ref": "green",
            "Allel1": "red",
            "Allel2": "orange",
            "HGSVC": "blue"
        }

        figure.for_each_trace(lambda trace: trace.update(marker=dict(color=color_mapping.get(trace.name, 'gray'))))

        unique_legend_names = set()
        for trace in figure.data:
            if trace.name in unique_legend_names:
                trace.showlegend = False
            elif trace.name in color_mapping:
                unique_legend_names.add(trace.name)
                trace.showlegend = True
            else:
                trace.showlegend = False

        name = "HGSVC" if st.session_state.analysis_mode == "indivisual sample" else st.session_state.analysis_mode 
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color='gray', size=20), name=name))

        figure.update_layout(
            title="Motif Occurrences",
            xaxis_title="Motif",
            yaxis_title="HGSVC vs Alleles",
        )

        xaxis_colors = {motif: motif_colors[idx] for idx, motif in enumerate(motif_names)}
        figure.update_xaxes(tickmode='array', tickvals=list(xaxis_colors.keys()), ticktext=[
            f'<span style="color:{xaxis_colors[motif]}">{motif}</span>' for motif in xaxis_colors.keys()
        ], tickangle=45)

        figure.update_yaxes(range=[0, df['Sample'].value_counts().max()])
        
        st.plotly_chart(figure, use_container_width=True)
    def plot_heatmap(self, combined_data, sort_by="Value"):
        

        if sort_by == "Value":
            x_sort = alt.EncodingSortField(field='Count', op='sum', order='descending')
            y_sort = alt.EncodingSortField(field='Count', op='sum', order='descending')
        else:
            x_sort = alt.SortField(field='Sample', order='ascending')
            y_sort = alt.SortField(field='Motif', order='ascending')
        combined_data = combined_data[combined_data['Motif'] != 'Interruption']
        heatmap = alt.Chart(combined_data).mark_rect().encode(
            x=alt.X('Sample:N', title='Sample', sort=x_sort, axis=alt.Axis(labelFontWeight='bold', labelColor='black', titleFontWeight='bold', titleColor='black')),
            y=alt.Y('Motif:N', title='Motif', sort=y_sort, axis=alt.Axis(labelFontWeight='bold', labelColor='black', titleFontWeight='bold', titleColor='black')),
            color=alt.Color('Count:Q', scale=alt.Scale(scheme='reds'), title='Count'),
            tooltip=['Sample', 'Motif', 'Count']
        ).configure_axis(
            labelFontWeight='bold',
            labelColor='black',
            titleFontWeight='bold',
            titleColor='black'
        ).configure_legend(
            labelFontWeight='bold',
            labelColor='black',
            titleFontWeight='bold',
            titleColor='black'
        )
        # Make the legend bold and black
        heatmap = heatmap.configure_legend(
            labelFontWeight='bold',
            labelColor='black',
            titleFontWeight='bold',
            titleColor='black'
        )

        # Remove interruptions from the heatmap data
        heatmap = heatmap.configure_legend(orient='top')
        st.altair_chart(heatmap, use_container_width=True)

    def bar_plot_motif_count(self, df, region, sort_by="Value"):
        df = df[df['Motif'] != "interruption"]
        
        total_copy_number = df.groupby('Sample').size().reset_index(name='Total Copy Number')
        total_copy_number = total_copy_number.sort_values(by='Total Copy Number', ascending=False)

        if sort_by == "Value":
            x_sort = alt.EncodingSortField(field='Total Copy Number', op='sum', order='descending')
        else:
            x_sort = alt.SortField(field='Sample', order='ascending')

        unique_samples = list(total_copy_number['Sample'].apply(lambda x: x.rsplit('_', 1)[0]).unique())
        color_palette = px.colors.qualitative.Vivid + px.colors.qualitative.Safe + px.colors.qualitative.Dark24 + px.colors.qualitative.Prism 

        color_mapping = {sample: color_palette[i] for i, sample in enumerate(unique_samples)}
        color_mapping = {sample: color_mapping[sample.rsplit('_', 1)[0]] for sample in total_copy_number['Sample']}

        bar_chart = alt.Chart(total_copy_number).mark_bar().encode(
            x=alt.X('Sample', sort=x_sort),
            y='Total Copy Number',
            tooltip=['Sample', 'Total Copy Number'],
            color=alt.Color('Sample', scale=alt.Scale(domain=list(color_mapping.keys()), range=list(color_mapping.values())))
        ).properties(
            width=600,
            height=400,
            title="Total Copy Number per Sample"
        ).configure_axis(
            labelFontWeight='bold',
            labelColor='black',
            titleFontWeight='bold',
            titleColor='black'
        ).configure_legend(
            labelFontWeight='bold',
            labelColor='black',
            titleFontWeight='bold',
            titleColor='black'
        )

        # Add the pathogenic threshold if applicable
        if region in st.session_state.pathogenic_TRs['region'].unique():
            pathogenic_threshold = st.session_state.pathogenic_TRs.loc[
                (st.session_state.pathogenic_TRs['region'] == region)
            ]['pathogenic_min'].values[0]

            # Create the threshold line
            threshold_line = alt.Chart(pd.DataFrame({'Total Copy Number': [pathogenic_threshold]})).mark_rule(
                color='red', strokeDash=[5, 5], size=3
            ).encode(
                y='Total Copy Number:Q'
            )

            threshold_pointer = alt.Chart(pd.DataFrame({'Total Copy Number': [pathogenic_threshold]})).mark_text(
                text='Pathogenic Threshold', align='left', dx=5, dy=-10, fontSize=12, color='black'
            ).encode(
                y='Total Copy Number:Q'
            )
            
            # Combine the bar chart, threshold line, and pointer
            bar_chart = alt.layer(bar_chart, threshold_line, threshold_pointer)

        # Apply configuration to the final chart
        bar_chart = bar_chart.configure_axisX(labelAngle=45, labelOverlap=False).configure_legend(orient='none', disable=True)

        # Render the chart
        st.altair_chart(bar_chart, use_container_width=True)

    
    def create_motif_dataframe(self,sequences, motif_colors, motif_ids, spans_list, motif_names):
        data = []
        interruptions_dict = set()
        for idx, sequence in enumerate(sequences):
            if sequence['sequence'] == "." or sequence['sequence'] == "":
                continue
            sequence_name = sequence['name']
            motif_ids_seq = motif_ids[idx]
            spans = spans_list[idx]
            ranges = self.parse_motif_range(spans)
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
        
        # print the interruptions as inturrption seen : 
        if interruptions_dict:
        
            interruptions_list = " | ".join([f"`{seq}`" for seq in interruptions_dict])
            st.markdown(f"**Interruptions Observed:** {interruptions_list}")
        else:
            st.subheader("No Significant Interruptions Detected")
        return pd.DataFrame(data)

    def get_color_palette(self,n):
        cmap = plt.get_cmap('tab20')  
        colors = [cmap(i) for i in range(n)]
        return ['#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in colors]
    def stack_plot(self, record, motif_names, sequences, span_list, motif_ids_list, sort_by="Value"):
        motif_colors = self.get_color_palette(len(record['motifs']))
        motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
        
        region = record['chr'] + ":" + str(record['pos'] - 1) + "-" + str(record['stop'] - 1)

        df = self.create_motif_dataframe(sequences, motif_colors, motif_ids_list, span_list, motif_names)
        if df.empty:
            return motif_colors, df
        df['Length'] = df['End'] - df['Start']
        
        default_height = 1500 
        chart_height = max(default_height, len(sequences) * 10)

        df['Sample'] = df['Sample'].apply(lambda x: x.replace("_pathogenic", ""))
        df['Sample'] = pd.Categorical(df['Sample'], categories=sorted(df['Sample'].unique()), ordered=True)
        df = df.sort_values(by=['Sample', 'Start', 'End'])
        df['Order'] = range(len(df))

        df_filtered = df[df['Motif'] != 'Interruption']

        for sample in df_filtered['Sample'].unique():
            df.loc[df['Sample'] == sample, 'Total Copy Number'] = df[df['Sample'] == sample].shape[0]

        min_copy_number = df['Total Copy Number'].min()
        max_copy_number = df['Total Copy Number'].max()

        pathogenic_threshold = 0
        gene_name = None
        inheritance = None
        disease = None
        if region in st.session_state.pathogenic_TRs['region'].unique():
            pathogenic_threshold = st.session_state.pathogenic_TRs.loc[
                (st.session_state.pathogenic_TRs['region'] == region)
            ]['pathogenic_min'].values[0]
            if pathogenic_threshold is None:
                pathogenic_threshold = 0
            gene_name = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == region
            ]['gene'].values[0]
            inheritance = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == region
            ]['inheritance'].values[0]
            disease = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == region
            ]['disease'].values[0]

        st.markdown(f"""
            <div style="display: flex; justify-content: space-between; font-size: 20px; color: #FF5733;">
                <div>
                    <strong>Minimum Copy Number:</strong> {min_copy_number}
                </div>
                <div>
                    <strong>Maximum Copy Number:</strong> {max_copy_number}
                </div>
            </div>
        """, unsafe_allow_html=True)

        if sort_by == "Value":
            y_sort = alt.EncodingSortField(field='Length', op='sum', order='descending')
        else:
            y_sort = alt.SortField(field='Sample', order='ascending')
        
        df["pathogenic"] = df["Total Copy Number"].apply(lambda x: "Pathogenic" if x >= pathogenic_threshold else "Not Pathogenic")
        pathogenic_threshold_length = pathogenic_threshold * len(motif_names[0])
        df['Sequence_length'] = df.groupby('Sample')['Length'].transform('sum')

        chart = alt.Chart(df).mark_bar().encode(
            y=alt.Y(
                'Sample', 
                sort=y_sort, 
                axis=alt.Axis(labelOverlap=False, ticks=False, labelColor='black', labelFontSize=12, labelFontWeight='bold', titleColor='black')  
            ),
            x=alt.X('Length', title='Length', stack='zero', axis=alt.Axis(labelColor='black', labelFontSize=12, labelFontWeight='bold', titleColor='black'), 
                    scale=alt.Scale(domain=[0, df['Sequence_length'].max() + 10])),
            color=alt.Color('Motif', scale=alt.Scale(domain=list(motif_names) + ['Interruption'], range=list(motif_colors.values()) + ['#FF0000'])),
            order=alt.Order('Order', sort='ascending'),
            tooltip=['Sample', 'Motif', 'Start', 'End', 'Sequence', 'pathogenic', 'Length', 'Sequence_length']
        ).properties(
            width=800,
            height=chart_height,
            title=alt.TitleParams(
                text="Motif occurrences across samples",
                anchor='middle',
                fontSize=20 
            )
        )

        # Adjust axis font size if chart height exceeds default
        if chart_height > default_height:
            chart = chart.configure_axisY(labelFontSize=10, titleFontSize=12)

        # Add gene info if available
        gene_info = []
        if gene_name:
            gene_info.append(f"Gene: {gene_name}")
        if inheritance:
            gene_info.append(f"Inheritance: {inheritance}")
        if disease:
            gene_info.append(f"Disease: {disease}")

        # Add pathogenic threshold rule if applicable
        if pathogenic_threshold > 0 and df.groupby('Sample')['Length'].sum().max() > pathogenic_threshold_length:
            rule = alt.Chart(pd.DataFrame({'x': [pathogenic_threshold_length], 'label': ['Pathogenic Threshold']})).mark_rule(
                color='red', strokeDash=[5, 5], clip=True
            ).encode(
                x=alt.X('x:Q', axis=alt.Axis(title='Length', labelColor='black', labelFontSize=12, labelFontWeight='bold', titleColor='black')),
                tooltip=['label', 'x'],
            ).encode(size=alt.value(5))

            arrow = alt.Chart(pd.DataFrame({'x': [pathogenic_threshold_length], 'y': [0], 'y2': [chart_height]})).mark_text(
                text='↑', align='left', baseline='bottom', dx=5, dy=45, fontSize=30, color='black', angle=290, clip=True
            ).encode(x=alt.X('x:Q'))

            arrow_line = alt.Chart(pd.DataFrame({'x': [pathogenic_threshold_length], 'y': [0], 'y2': [chart_height]})).mark_rule(
                color='black', clip=True
            ).encode(x=alt.X('x:Q'))

            arrow_text = alt.Chart(pd.DataFrame({'x': [pathogenic_threshold_length], 'y': [chart_height / 2], 'text': ['Pathogenic Threshold']})).mark_text(
                align='left', baseline='middle', dx=10, dy=25, fontSize=12, color='black'
            ).encode(x='x:Q', text='text')

            combined_chart = alt.layer(chart, rule, arrow, arrow_line, arrow_text).resolve_scale(x='shared')
            
            combined_chart = combined_chart.properties(
                title=alt.TitleParams(
                    text="Motif occurrences across samples",
                    subtitle=gene_info,
                    anchor='middle',
                    fontSize=20,
                    subtitleFontSize=15,
                    subtitleColor='green',
                    subtitlePadding=20  
                )
            ).configure_axis(
                labelFontSize=10,
                titleFontSize=12,
            ).configure_legend(
                labelFontSize=12,
                labelFontWeight='bold',
                labelColor='black',
                titleFontSize=14,
                titleFontWeight='bold',
                titleColor='black'
            ).properties(
                padding={'left': 10, 'right': 50, 'top': 30, 'bottom': 10}
            )

            st.altair_chart(combined_chart, use_container_width=True)
        else:
            chart = chart.properties(
                title=alt.TitleParams(
                    text="Motif occurrences across samples",
                    subtitle=gene_info,
                    anchor='middle',
                    fontSize=20,
                    subtitleFontSize=15,
                    subtitleColor='green',
                    subtitlePadding=20,
                )
            ).configure_legend(
                labelFontSize=12,
                labelFontWeight='bold',
                labelColor='black',
                titleFontSize=14,
                titleFontWeight='bold',
                titleColor='black'
            )
            st.altair_chart(chart, use_container_width=True)

        return motif_colors, df

    def count_motifs(self, motif_ids):
        motif_count = {}
        for idx, motif in enumerate(motif_ids):
            if motif in motif_count:
                motif_count[motif] += 1
            else:
                motif_count[motif] = 1
        
        return motif_count



    def visulize_TR_with_dynamic_sequence(self,record,hgsvc_records, left_column, right_column,motif_colors,CN1_col,CN2_col, show_comparison):
        
        if record['motif_ids_h1'] == ['.'] and record['motif_ids_h2'] == ['.']:
            st.info("No motifs found in the region")
            return
        motif_names = record['motifs']
        reference_copy_number = record['ref_CN']
        motif_count_ref = self.count_motifs(record['motif_ids_ref'])
        found_motifs_ref = list(motif_count_ref.keys())
        found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]

        motif_count_h1 = self.count_motifs(record['motif_ids_h1'])
        found_motifs_h1 = list(motif_count_h1.keys())
        found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
        motif_count_h2 = self.count_motifs(record['motif_ids_h2'])
        found_motifs_h2 = list(motif_count_h2.keys())
        found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
        motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
        motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
        total_copy_number_h1 = str(record['spans'][1]).count('-')
        CN1_col.markdown(f"""
                <div style="font-size: 20px; color: #FF5733;">
                    <strong>Allele 1 Total copy number:</strong> {total_copy_number_h1}
                </div>
            """, unsafe_allow_html=True)
        if record['alt_allele2'] != '':
                total_copy_number_h2 = str(record['spans'][2]).count('-')
                CN2_col.markdown(f"""
                    <div style="font-size: 20px; color: #FF5733;">
                        <strong>Allele 2 Total copy number:</strong> {total_copy_number_h2}
                    </div>
                """, unsafe_allow_html=True)   
        with right_column:
            self.display_motif_legend(motif_names, motif_colors, right_column)
        with left_column:
            tab1, tab2,tab3 = st.tabs(["Alleles", "Alleles vs Ref", "Alleles vs Pop"])
            with tab2:
                self.display_dynamic_sequence_with_highlighted_motifs("Ref", record['ref_allele'], record['motif_ids_ref'], record['spans'][0], motif_colors, motif_names)
                alt_allele1 = record['alt_allele1']
                self.display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names)
                if record['alt_allele2'] != '':
                    alt_allele2 = record['alt_allele2'] 
                    self.display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names)
            with tab1:
                alt_allele1 = record['alt_allele1']
                self.display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names)

                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        plot_motif_bar(motif_count_h1, motif_names, motif_colors)

                if record['alt_allele2'] != '':

                    alt_allele2 = record['alt_allele2'] 
                    self.display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names)
                    
                    plot_container_h2 = st.empty()
                    with plot_container_h2:
                        if show_comparison == False:
                            plot_motif_bar(motif_count_h2, motif_names, motif_colors)
            with tab3:
                if hgsvc_records:
                    self.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
                else:
                    st.info("no population data found")


    def display_motif_legend(self, motifs, motif_colors, right_column):
        st.markdown("### Motif Legend")
        st.markdown('<div style="max-height:400px; overflow-y:scroll;">', unsafe_allow_html=True)  
        if  isinstance(motifs, tuple):
            motifs = list(motifs)
        elif not isinstance(motifs, list):
            motifs = [motifs]
        for idx, motif in enumerate(motifs):
            color = motif_colors[idx]
            st.markdown(
                f'<div id="legend-motif-{idx}" class="legend-item motif-{idx}" style="background-color:{color};color:white;padding:5px;margin-bottom:10px;border-radius:5px;'
                f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);'
                f' white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="{motif}">'
                f' Motif {idx}: {motif}</div>', unsafe_allow_html=True)
        st.markdown(
            f'<div id="legend-motif-interruption" class="legend-item motif-interruption" style="background-color:#FF0000;color:black;padding:5px;margin-bottom:10px;border-radius:5px;'
            f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
            f' Interruption</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)







# --- Helper Functions Below ---

# def get_simple_spans(spans_string):
#     """
#     Take a string like '0(0-27)_1(0-28)' and pull out just the span coordinates as strings,
#     e.g. ['(0-27)', '(0-28)']. Ignores motif ids.
#     """
#     if not spans_string:
#         return []
#     simple_spans = []
#     for motif_span in spans_string.split('_'):
#         match = re.match(r'(?:\d+)?\((\d+-\d+)\)', motif_span)
#         if match:
#             span = match.group(1)
#             simple_spans.append(f'({span})')
#     return simple_spans

# def split_spans_by_occurrences(spans_string, motif_ids, motifs):
#     """
#     Split all TRGT-style (multi-span) regions into subspans per-motif occurrence for visualization.
#     Handles multiple spans and motif counts, e.g.:
#     spans_string: "15(0-51)_0(51-102)_9(102-155)_..."
#     motif_ids: [1, 1, 1, ...] (total n)
#     Returns a flattened list of per-motif subspans, e.g. ['(0-51)', '(51-102)', ...]
#     """
#     import re

#     if not spans_string or not motif_ids or not isinstance(motif_ids, list):
#         return []

#     # Parse all spans in the TRGT string
#     # Format: "motif_idx(start-end)_motif_idx(start-end)_..."
#     span_matches = re.findall(r'(?:\d+)?\((\d+)-(\d+)\)', spans_string)
#     if not span_matches or len(motif_ids) == 0:
#         return []

#     subspans = []
#     for idx, span in enumerate(span_matches):

#         span_start = int(span[0])
#         span_end = int(span[1])
#         raw_len = span_end - span_start 
#          # Always run at least once, but can be 0
#         # Handle unit length: motif spanning single bases
#         unit_len = len(motifs[motif_ids[idx]])
#         actual_count = max(1, raw_len // unit_len)
#         curr_start = span_start
#         for repeat in range(actual_count):
#             curr_end = curr_start + unit_len
#             subspans.append(f'({curr_start}-{curr_end})')
#             curr_start = curr_end
#     subspans = ''.join(subspans)
#     return subspans
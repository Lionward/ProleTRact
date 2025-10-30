import streamlit as st
import altair as alt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import re 
import pysam
import numpy as np
import importlib
display_dynamic_sequence_with_highlighted_motifs = importlib.reload(__import__("vis_helper")).display_dynamic_sequence_with_highlighted_motifs
display_motifs_as_bars = importlib.reload(__import__("vis_helper")).display_motifs_as_bars
motif_legend_html = importlib.reload(__import__("vis_helper")).motif_legend_html
plot_motif_bar = importlib.reload(__import__("vis_helper")).plot_motif_bar
interpret_genotype = importlib.reload(__import__("vis_helper")).interpret_genotype
display_genotype_card = importlib.reload(__import__("vis_helper")).display_genotype_card
display_genotype_badge = importlib.reload(__import__("vis_helper")).display_genotype_badge
create_genotype_comparison_matrix = importlib.reload(__import__("vis_helper")).create_genotype_comparison_matrix


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
            if st.session_state.cohort_mode == "assembly":
                record = self.parse_record_assembly(files[i], region)
            else:
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
                'gt': '',
                'id': '',
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
            'gt': str(rec.samples[0]['GT'][0]),
            'id': rec.id,
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
    def render_region_display(self, markdown_placeholder, region):
        html = f"""
            <style>
                .region-badge {{
                    display: inline-block;
                    background: #ECEAFB;
                    padding: 6px 12px;
                    border-radius: 50px;
                    font-size: 8px;
                    font-weight: 700;
                    margin-bottom: 20px;
                    backdrop-filter: blur(10px);
                    border: 1px solid #C9BEEF;
                    letter-spacing: 1px;
                }}
            </style>
            
            <div style="display: flex; justify-content: center; align-items: center; min-height: 80px;">
                <div style="
                    display: flex; 
                    align-items: center; 
                    background: #ECEAFB; 
                    padding: 9px 16px; 
                    border-radius: 20px; 
                    box-shadow: 0 4px 16px 2px rgba(118, 75, 162, 0.12); 
                    font-size: 1.5rem; 
                    font-weight: 700;
                    border: 1px solid #C9BEEF;
                    ">
                    <span style="width: 20px; height: 20px; background: #764ba2; border-radius: 50%; margin-right: 20px; box-shadow: 0 0 8px 2px #A184D6;"></span>
                    Region: {region}
                </div>
            </div>
        """
        markdown_placeholder.html(html)

    def visulize_cohort(self):
        if 'cohorts_records_map' in st.session_state:
            region_options = list(st.session_state.cohorts_records_map.values())
            
            # Safe index access with bounds checking
            regions_idx = st.session_state.get('regions_idx', 0)
            if regions_idx >= len(region_options):
                regions_idx = 0
                st.session_state.regions_idx = 0
            
            default_region = region_options[regions_idx] if region_options else ""
            
            # Cache the full options list to avoid recreating it
            if 'cached_region_options_cohort' not in st.session_state:
                st.session_state.cached_region_options_cohort = region_options
            
            # Single unified field: text input with autocomplete suggestions
            # Track if a selection was made
            if 'region_selected_cohort' not in st.session_state:
                st.session_state.region_selected_cohort = ""
            
            search_query = st.sidebar.text_input(
                "🔍 Search region:", 
                value=st.session_state.region_selected_cohort,
                key="region_search_cohort",
                help="Type to search and filter results",
                placeholder="Type to search..."
            )
            
            # Use markdown+CSS trick to make selectbox text black
            st.markdown("""
                <style>
                /* Ensure all selectbox texts are black, regardless of state */
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] span {
                    color: black !important;
                }
                [data-testid="stSidebar"] .stSelectbox label, 
                [data-testid="stSidebar"] .stSelectbox div[role="listbox"] span {
                    color: black !important;
                }
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] input,
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="combobox"] span,
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="button"] span,
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="option"] span,
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="listbox"] span {
                    color: black !important;
                }
                /* Also set all the selectbox selected value text to black */
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] input {
                    color: black !important;
                    font-weight: 700;
                }
                /* For v1.31+ (possible dark mode or material theme changes): */
                [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[aria-selected="true"] span {
                    color: black !important;
                }
                </style>
            """, unsafe_allow_html=True)
    
            
            # Filter and show suggestions only when typing
            if search_query:
                search_lower = search_query.lower()
                filtered = [r for r in st.session_state.cached_region_options_cohort if search_lower in r.lower()]
                filtered_regions = filtered[:10]  # Limit to 10 results
                total_matches = len(filtered)
                
                # Show filtered suggestions as a dropdown without label to feel unified
                if filtered_regions:
                    region = st.sidebar.selectbox(
                        " ",  # Empty label
                        filtered_regions, 
                        index=0,
                        key="region_suggest",
                        help="Select from suggestions",
                        label_visibility="collapsed"
                    )
                    # Update the text input with the selected region
                    st.session_state.region_selected_cohort = region
                else:
                    region = search_query
                    
                # Show match count
                if total_matches > 10:
                    st.sidebar.markdown(f"<span style='font-size:11px; color:orange;'>Showing 10 of {total_matches:,} matches</span>", unsafe_allow_html=True)
                else:
                    st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>{total_matches:,} matches</span>", unsafe_allow_html=True)
            else:
                # When not typing, use the current region from session state
                region = default_region
                total_matches = len(st.session_state.cached_region_options_cohort)
                st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>Search to find from {total_matches:,} regions</span>", unsafe_allow_html=True)

            if 'regions_idx' not in st.session_state:
                st.session_state.regions_idx = 0
            
            # Beautiful navigation buttons in sidebar
            st.sidebar.markdown("""
                <div style="
                    display: flex;
                    align-items: center;
                    gap: 10px;
                    background: linear-gradient(92deg, #667eea 0%, #a084e8 100%);
                    padding: 10px 18px;
                    border-radius: 14px;
                    box-shadow: 0 2px 14px rgba(102, 126, 234, 0.10);
                    margin-bottom: 8px;
                    ">
                    <span style="
                        font-size: 28px;
                        color: #4338ca;
                        margin-right: 6px;
                        filter: drop-shadow(0 2px 6px rgba(65,0,140,0.18));
                        ">🧭</span>
                    <span style="
                        font-size: 1.25rem; 
                        font-weight: 700; 
                        letter-spacing: 0.03em; 
                        color: #312e81;
                        ">Region Navigation</span>
                </div>
            """, unsafe_allow_html=True)
            nav_col1, nav_col2 = st.sidebar.columns(2, gap="small")
            with nav_col1:
                if st.button("◀ Previous", use_container_width=True, key="prev_region"):
                    region = None
                    st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)
            with nav_col2:
                if st.button("Next ▶", use_container_width=True, key="next_region"):
                    region = None
                    st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len( st.session_state.cohorts_records_map )-1)
            
            # Add beautiful styling for sidebar navigation buttons
            st.markdown("""
                <style>
                    /* Beautiful sidebar navigation buttons */
                    [data-testid="stSidebar"] button[kind="secondary"] {
                        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                        color: white !important;
                        border: none !important;
                        border-radius: 12px !important;
                        font-weight: 700 !important;
                        font-size: 15px !important;
                        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
                        box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3) !important;
                    }
                    [data-testid="stSidebar"] button[kind="secondary"]:hover {
                        transform: translateY(-2px) scale(1.05) !important;
                        box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4) !important;
                        background: linear-gradient(135deg, #764ba2 0%, #667eea 100%) !important;
                    }
                </style>
                <script>
                (function() {
                    function styleNavButtons() {
                        const sidebarButtons = document.querySelectorAll('[data-testid="stSidebar"] button');
                        sidebarButtons.forEach(btn => {
                            const text = btn.textContent || btn.innerText || '';
                            if (text.includes('◀ Previous') || text.includes('Next ▶')) {
                                btn.style.background = 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)';
                                btn.style.color = 'white';
                                btn.style.border = 'none';
                                btn.style.borderRadius = '12px';
                                btn.style.fontWeight = '700';
                                btn.style.fontSize = '15px';
                                btn.style.transition = 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)';
                                btn.style.boxShadow = '0 4px 15px rgba(102, 126, 234, 0.3)';
                            }
                        });
                    }
                    styleNavButtons();
                    setTimeout(styleNavButtons, 100);
                    setTimeout(styleNavButtons, 500);
                    const observer = new MutationObserver(styleNavButtons);
                    observer.observe(document.body, { childList: true, subtree: true });
                })();
                </script>
            """, unsafe_allow_html=True)
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
                #st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len( st.session_state.cohorts_records_map )-1)
                region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
            mode_placeholder = st.empty()
            
            # Create container for region display
            middel, _ = st.columns([1, 0.1], gap="small")
            self.render_region_display(middel, region)
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

        gt = rec.samples[0]['GT']
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
            'motif_ids_ref': rec.info['MOTIF_IDs_REF'].split("_"),
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
        # For debugging/inspection purposes
        return record


    def get_color_palette(self, num_colors):
        """
        Generate a color palette optimized for motif visualization with good distinctiveness
        """
        # Primary palette matching your app's gradient theme
        primary_palette = [
            '#667eea',  # Primary blue
            '#764ba2',  # Primary purple
            '#4fd1c7',  # Teal
            '#68d391',  # Green
            '#f6ad55',  # Orange
            '#fc8181',  # Coral
            '#f093fb',  # Pink
            '#f6e05e',  # Yellow
        ]
        
        # Secondary palette - darker shades
        secondary_palette = [
            '#5a67d8', '#6b46c1', '#38a169', '#dd6b20',
            '#e53e3e', '#d53f8c', '#d69e2e', '#319795'
        ]
        
        # Tertiary palette - lighter shades
        tertiary_palette = [
            '#a3bffa', '#b794f4', '#9ae6b4', '#fbd38d',
            '#fed7d7', '#fbb6ce', '#faf089', '#81e6d9'
        ]
        
        # Combine all palettes
        combined_palette = primary_palette + secondary_palette + tertiary_palette
        
        # If we need more colors, generate them algorithmically
        if num_colors > len(combined_palette):
            import colorsys
            import math
            
            generated_colors = []
            # Use golden ratio for optimal color distribution
            golden_ratio_conjugate = (math.sqrt(5) - 1) / 2
            
            for i in range(num_colors - len(combined_palette)):
                # Generate colors with good separation in HSV space
                hue = (i * golden_ratio_conjugate) % 1.0
                # Vary saturation and value for better distinction
                saturation = 0.6 + 0.2 * (i % 3) / 3
                value = 0.7 + 0.2 * ((i + 1) % 2)
                
                rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                hex_color = '#{:02x}{:02x}{:02x}'.format(
                    int(rgb[0] * 255), 
                    int(rgb[1] * 255), 
                    int(rgb[2] * 255)
                )
                generated_colors.append(hex_color)
            
            return combined_palette + generated_colors[:num_colors]
        
        return combined_palette[:num_colors]
        
      
    def visulize_region(self):
        if 'regions_idx' not in st.session_state:
            st.session_state.regions_idx = 0
        region_options = list(st.session_state.records_map.values())
        
        # Safe index access with bounds checking
        regions_idx = st.session_state.get('regions_idx', 0)
        if regions_idx >= len(region_options):
            regions_idx = 0
            st.session_state.regions_idx = 0
        
        default_region = region_options[regions_idx] if region_options else ""
        st.sidebar.markdown("### Select Region to Visualize")
        
        # Cache the full options list to avoid recreating it
        if 'cached_region_options' not in st.session_state:
            st.session_state.cached_region_options = region_options
        
        # Track if a selection was made
        if 'region_selected_ind' not in st.session_state:
            st.session_state.region_selected_ind = ""
        
        # Searchable field: text input that filters selectbox
        search_query = st.sidebar.text_input(
            "🔍 Search region:", 
            value=st.session_state.region_selected_ind,
            key="region_search",
            help="Type to search and filter results",
            placeholder="Type to search..."
        )
        
        # Filter options based on search and show only 10 results
        if search_query:
            search_lower = search_query.lower()
            filtered = [r for r in st.session_state.cached_region_options if search_lower in r.lower()]
            filtered_regions = filtered[:10]  # Show only 10 results
            total_matches = len(filtered)
            
            # Show the filtered selectbox without label to feel unified
            if filtered_regions:
                region = st.sidebar.selectbox(
                    " ",  # Empty label
                    filtered_regions, 
                    index=0,
                    key="region_suggest",
                    help="Select from filtered results",
                    label_visibility="collapsed"
                )
                # Update the text input with the selected region
                st.session_state.region_selected_ind = region
            else:
                region = search_query
                
            if total_matches > 10:
                st.sidebar.markdown(f"<span style='font-size:11px; color:orange;'>Showing 10 of {total_matches:,} matches</span>", unsafe_allow_html=True)
            else:
                st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>{total_matches:,} matches</span>", unsafe_allow_html=True)
        else:
            # No search - show default region
            region = default_region
            total_matches = len(st.session_state.cached_region_options)
            st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>Search to find from {total_matches:,} regions</span>", unsafe_allow_html=True)
        
        display_option = st.sidebar.radio("Select Display Type", 
                                            ("Sequence with Highlighted Motifs", "Bars"))

        # Beautiful navigation buttons in sidebar
        st.sidebar.markdown("### 🧭 Navigation")
        nav_col1, nav_col2 = st.sidebar.columns(2, gap="small")
        with nav_col1:
            if st.button("◀ Previous", use_container_width=True, key="prev_individual"):
                region = None
                st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

        with nav_col2:
            if st.button("Next ▶", use_container_width=True, key="next_individual"):
                region = None
                st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.records_map) - 1)
        
        # Add beautiful styling for sidebar navigation buttons
        st.markdown("""
            <style>
                /* Beautiful sidebar navigation buttons */
                [data-testid="stSidebar"] button[kind="secondary"] {
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                    color: white !important;
                    border: none !important;
                    border-radius: 12px !important;
                    font-weight: 700 !important;
                    font-size: 15px !important;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
                    box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3) !important;
                }
                [data-testid="stSidebar"] button[kind="secondary"]:hover {
                    transform: translateY(-2px) scale(1.05) !important;
                    box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4) !important;
                    background: linear-gradient(135deg, #764ba2 0%, #667eea 100%) !important;
                }
            </style>
            <script>
            (function() {
                function styleNavButtons() {
                    const sidebarButtons = document.querySelectorAll('[data-testid="stSidebar"] button');
                    sidebarButtons.forEach(btn => {
                        const text = btn.textContent || btn.innerText || '';
                        if (text.includes('◀ Previous') || text.includes('Next ▶')) {
                            btn.style.background = 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)';
                            btn.style.color = 'white';
                            btn.style.border = 'none';
                            btn.style.borderRadius = '12px';
                            btn.style.fontWeight = '700';
                            btn.style.fontSize = '15px';
                            btn.style.transition = 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)';
                            btn.style.boxShadow = '0 4px 15px rgba(102, 126, 234, 0.3)';
                        }
                    });
                }
                styleNavButtons();
                setTimeout(styleNavButtons, 100);
                setTimeout(styleNavButtons, 500);
                const observer = new MutationObserver(styleNavButtons);
                observer.observe(document.body, { childList: true, subtree: true });
            })();
            </script>
        """, unsafe_allow_html=True)
        
        middel, spacer = st.columns([1, 0.1], gap="small")
        REF, CN1_col, CN2_col = st.columns([1, 1, 1])
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

        self.render_region_display(middel, st.session_state.records_map[st.session_state.regions_idx])


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

        container = st.container()
        motif_colors = self.get_color_palette(len(record['motifs']))
        motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
            
        col1,col2 = st.sidebar.columns([1,1])



        if display_option == "Sequence with Highlighted Motifs":
            self.visulize_TR_with_dynamic_sequence(record,hgsvc_records, container,motif_colors,CN1_col,CN2_col, st.session_state.get('show_comparison', False))

        elif display_option == "Bars":
            self.display_motifs_with_bars(record, container ,motif_colors,CN1_col,CN2_col, st.session_state.get('show_comparison', False),hgsvc_records)


    def parse_motif_range(self,motif_range):
        pattern = re.compile(r'\((\d+)-(\d+)\)')
        matches = pattern.findall(motif_range)
        ranges = [(int(start)-1, int(end)-1) for start, end in matches]
        return ranges
    
    
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




        
    def display_motifs_with_bars(self, record, containter, motif_colors, CN1_col, CN2_col, show_comparison, hgsvc_records):
        motif_names, motif_count_h1, motif_count_h2 = self.parse_motif_in_region(record)
        if motif_count_h1 == None and motif_count_h2 == None:
            st.info("No motifs found in the region")
            return
        st.markdown("""
            <style>
                .stTabs [data-baseweb="tab-list"] {
                    gap: 4px;
                    background: #f8fafc;
                    padding: 2px 2px;
                    border-radius: 12px;
                    margin-bottom: 12px;
                    min-height: 15px;
                }
                .stTabs [data-baseweb="tab"] {
                    height: 18px;
                    min-height: 18px;
                    min-width: 28px;
                    white-space: pre-wrap;
                    background: white;
                    border-radius: 8px;
                    border: 1px solid #e2e8f0;
                    color: #64748b;
                    font-weight: 600;
                    font-size: 10px;
                    padding: 0 7px;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                    line-height: 1.1;
                }
                .stTabs [data-baseweb="tab"]:hover {
                    background: #f1f5f9;
                    border-color: #cbd5e1;
                    color: #475569;
                    transform: translateY(-1.6px);
                    box-shadow: 0 2px 8px rgba(0,0,0,0.10);
                }
                .stTabs [aria-selected="true"] {
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                    color: white !important;
                    border-color: #667eea !important;
                    box-shadow: 0 4px 14px rgba(102, 126, 234, 0.2) !important;
                    transform: translateY(-1.6px);
                }
                .stTabs [data-baseweb="tab-highlight"] {
                    background-color: transparent !important;
                }
            </style>
        """, unsafe_allow_html=True)

        

        with containter:
            tab1, tab2, tab3 = st.tabs([
                "🧬 **Alleles**", 
                "🔄 **Alleles vs Ref**", 
                "🌍 **Alleles vs Pop**"
            ])
            st.markdown(
                "<style>.tab-content {font-size: 20px;}</style>",
                unsafe_allow_html=True,
            )


            with tab2:
                motif_legend_html(record['motif_ids_ref'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("### Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                display_motifs_as_bars("Ref", motif_colors, record['motif_ids_ref'], record['spans'][0], record['ref_allele'], motif_names, None)
                display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names, record['supported_reads_h1'])
                if record['alt_allele2'] != '':
                    display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names, record['supported_reads_h2'])
            with tab1:
                motif_legend_html(record['motif_ids_h1'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("###  Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names, record['supported_reads_h1'])
                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        with st.expander("Show motif bar plot for Allel1", expanded=False):
                            plot_motif_bar(motif_count_h1, motif_names, motif_colors)
                
                if record['alt_allele2'] != '':
                    display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names, record['supported_reads_h2'])
                    plot_container_h2 = st.empty()

                    with plot_container_h2:
                        if show_comparison == False:
                            with st.expander("Show motif bar plot for Allel2", expanded=False):
                                plot_motif_bar(motif_count_h2, motif_names, motif_colors)

            with tab3:
                if hgsvc_records:
                    # Add genotype comparison for population data
                    st.markdown("### 🧬 Population Genotype Comparison")
                    population_genotypes = {"Current Sample": record['gt']}
                    for sample_name, sample_record in hgsvc_records.items():
                        if 'gt' in sample_record:
                            population_genotypes[sample_name] = sample_record['gt']
                    
                    if len(population_genotypes) > 1:
                        create_genotype_comparison_matrix(population_genotypes)
                        st.markdown("---")
                    
                    self.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
                else:
                    st.info("no population data found")

            



    def plot_Cohort_results(self,cohort_records):
        sequences = []
        span_list = []
        motif_ids_list = []
        # make space between the last print 
        sort_by = st.radio("Sort by:", ("Value", "Sample Name"), horizontal=True, key="sort_by_cohort")
        
        # Extract genotype information for comparison
        genotypes_dict = {}
        for key in cohort_records.keys():
            genotypes_dict[key] = cohort_records[key]['gt']
        
        # Display genotype comparison matrix
        st.markdown("---")
        create_genotype_comparison_matrix(genotypes_dict)
        st.markdown("---")
        if st.session_state.cohort_mode == "assembly":
            for key in cohort_records.keys():
                sequences.append({'name': f'{key}', 'sequence': cohort_records[key]['alt_allele']})
                span_list.append(cohort_records[key]['spans'])
                motif_ids_list.append(cohort_records[key]['motif_ids_h'])
        else:
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
        
        # Enhanced radio button styling
        # st.markdown("""
        #     <style>
        #         .stRadio > div {
        #             background: rgba(255, 255, 255, 0.9);
        #             padding: 15px;
        #             border-radius: 16px;
        #             border: 1px solid rgba(102, 126, 234, 0.2);
        #             margin-bottom: 20px;
        #         }
        #         .stRadio label {
        #             font-weight: 600 !important;
        #             color: #4a5568 !important;
        #         }
        #     </style>
        # """, unsafe_allow_html=True)
        
        # Enhanced radio group display using Streamlit columns for a visually appealing "Sort by"
        st.html("""
            <style>
                /* Outer radio styling */
                .custom-radio-group {
                    display: flex;
                    gap: 22px;
                    align-items: center;
                    background: rgba(255, 255, 255, 0.97);
                    border-radius: 16px;
                    padding: 18px 32px 18px 20px;
                    border: 1.5px solid #c7d2fe;
                    box-shadow: 0 2px 16px rgba(102, 126, 234, 0.07);
                    margin: 22px 0 18px 0;
                    backdrop-filter: blur(14px);
                }
                .custom-radio-label {
                    display: flex !important;
                    align-items: center;
                    gap: 7px;
                    margin-right: 0 !important;
                }
                .custom-radio-btn {
                    appearance: none;
                    width: 22px;
                    height: 22px;
                    border-radius: 50%;
                    border: 2.4px solid #bfc8fa;
                    transition: border .2s, box-shadow .2s;
                    outline: none;
                    margin-right: 12px;
                    box-sizing: border-box;
                    background: white;
                    box-shadow: 0 0 0px 4px #e1eaff00;
                    position: relative;
                }
                .custom-radio-btn:checked {
                    border-color: #667eea;
                    background: linear-gradient(135deg, #667eea 60%, #764ba2 100%);
                    box-shadow: 0 0 0px 4px #bfc8fa44;
                }
                .custom-radio-btn:checked::before {
                    content: '';
                    position: absolute;
                    left: 6.5px; top: 6.5px;
                    width: 7px; height: 7px;
                    border-radius: 50%;
                    background: #fff;
                }
                .custom-radio-label.active,
                .custom-radio-label:hover {
                    background: linear-gradient(135deg, #f3f5fa 90%, #eef0fe 100%);
                    border-radius: 9px;
                    box-shadow: 0 2px 12px #e1eaff33;
                    cursor: pointer;
                }
                .custom-radio-label .radio-emoji {
                    font-size: 1.16em;
                }
                .sortby-title {
                    font-size: 19px;
                    font-weight: 700;
                    color: #43486b;
                    margin-right: 9px;
                }
            </style>
        """)

        sortby_options = [
            {"label": '<span class="radio-emoji">📊</span> Value', "value": "📊 Value"},
            {"label": '<span class="radio-emoji">👥</span> Sample Name', "value": "👥 Sample Name"}
        ]

        # Maintain radio state with Streamlit session_state
        if "sort_by_cohort" not in st.session_state:
            st.session_state["sort_by_cohort"] = sortby_options[0]["value"]

        st.html(
            '<div class="custom-radio-group">' +
            '<span class="sortby-title">Sort by:</span>' +
            ''.join(
                f'''<label class="custom-radio-label{' active' if st.session_state["sort_by_cohort"] == option['value'] else ''}">
                    <input type="radio" class="custom-radio-btn" name="sort_by_cohort" value="{option['value']}" 
                        {'checked' if st.session_state["sort_by_cohort"] == option['value'] else ''} 
                        onclick="window.dispatchEvent(new CustomEvent('sortbyChange', {{detail: '{option['value']}'}}))">
                    {option['label']}
                </label>''' 
                for option in sortby_options
            ) +
            '</div>',
            
        )

        # Streamlit custom events solution to handle client-side radio change
        st.components.v1.html(f"""
        <script>
            window.addEventListener('sortbyChange', function(e) {{
                // Use a cookie to communicate - Streamlit will reload when widget is redrawn
                document.cookie = "sortby_cohort_choice=" + encodeURIComponent(e.detail) + "; path=/";
                window.location.reload();
            }});
        </script>
        """, height=0)

        motif_colors = self.get_color_palette(len(motif_names))

        # Server-side get cookie function
        def get_cookie_value(key):
            try:
                cookies_param = st.query_params().get('__streamlit_cookies__')
                if cookies_param:
                    # Handle both string and list cases
                    if isinstance(cookies_param, list):
                        cookies = cookies_param[0] if cookies_param else ''
                    else:
                        cookies = str(cookies_param)
                else:
                    cookies = ''
                m = _re.search(rf'{_re.escape(key)}=(.*?)(?:;|$)', cookies)
                return m.group(1) if m else None
            except Exception:
                return None

        # Streamlit doesn't provide cookie parsing directly, so we try with components
        sortby_cookie = get_cookie_value('sortby_HGSVC_choice')
        if sortby_cookie in [o["value"] for o in sortby_options]:
            st.session_state["sort_by_HGSVC"] = sortby_cookie

        sort_by = st.session_state.get("sort_by_HGSVC", sortby_options[0]["value"])
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
        
        # Create two columns for side-by-side plots


        st.markdown('<div class="plot-card">', unsafe_allow_html=True)
        self.bar_plot_motif_count(df, region, sort_by)
        st.markdown('</div>', unsafe_allow_html=True)
        
        
        # Scatter plot in full width below
        st.markdown('<div class="plot-card-full">', unsafe_allow_html=True)
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
                    marker=dict(
                        size=20,
                        line=dict(width=2, color='white')
                    ),
                    hovertemplate=f"<b>{sample}</b><br>Motif: %{{x}}<br>Count: %{{y}}<extra></extra>"
                )
                if sample in ["Ref", "Allel1", "Allel2"]:
                    if sample == "Ref":
                        trace['marker']['symbol'] = "diamond"
                        ref_data.append(trace)
                    else:
                        trace['marker']['symbol'] = "triangle-down"
                        allele_data.append(trace)
                else:
                    figure.add_trace(trace)

        for trace in ref_data + allele_data:
            figure.add_trace(trace)

        color_mapping = {
            "Ref": "#10B981",
            "Allel1": "#EF4444", 
            "Allel2": "#F59E0B",
            "HGSVC": "#3B82F6"
        }

        figure.for_each_trace(lambda trace: trace.update(
            marker=dict(
                color=color_mapping.get(trace.name, '#6B7280'),
                line=dict(width=2, color='white')
            )
        ))

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
        figure.add_trace(go.Scatter(
            x=[None], y=[None], 
            mode='markers', 
            marker=dict(color='#6B7280', size=20, line=dict(width=2, color='white')), 
            name=name
        ))

        figure.update_layout(
            title=dict(
                text="Motif Occurrences Across Samples",
                x=0.2,
                font=dict(size=20, color='#1F2937')
            ),
            xaxis_title=dict(text="Motif", font=dict(size=14, color='#4B5563')),
            yaxis_title=dict(text="Occurrence Count", font=dict(size=14, color='#4B5563')),
            plot_bgcolor='rgba(255,255,255,0.9)',
            paper_bgcolor='rgba(255,255,255,0.9)',
            hoverlabel=dict(
                bgcolor="white",
                font_size=12,
                font_family="Inter"
            ),
            legend=dict(
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='rgba(0,0,0,0.1)',
                borderwidth=1
            )
        )

        xaxis_colors = {motif: motif_colors[idx] for idx, motif in enumerate(motif_names)}
        figure.update_xaxes(
            tickmode='array', 
            tickvals=list(xaxis_colors.keys()), 
            ticktext=[
                f'<span style="color:{xaxis_colors[motif]}; font-weight:600">{motif}</span>' for motif in xaxis_colors.keys()
            ], 
            tickangle=45,
            gridcolor='rgba(0,0,0,0.1)'
        )

        figure.update_yaxes(
            range=[0, df['Sample'].value_counts().max()],
            gridcolor='rgba(0,0,0,0.1)'
        )
        
        st.plotly_chart(figure, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)


    def bar_plot_motif_count(self, df, region, sort_by="Value"):
        df = df[df['Motif'] != "interruption"]
        
        total_copy_number = df.groupby('Sample').size().reset_index(name='Total Copy Number')
        total_copy_number = total_copy_number.sort_values(by='Total Copy Number', ascending=False)
        if sort_by == "Value":
            x_sort = alt.EncodingSortField(field='Total Copy Number', op='sum', order='descending')
        else:
            x_sort = alt.SortField(field='Sample', order='ascending')

        unique_samples = list(total_copy_number['Sample'].apply(lambda x: x.rsplit('_', 1)[0]).unique())
        
        # Updated color palette to match the theme
        color_palette = [
            '#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391',
            '#f6e05e', '#f6ad55', '#fc8181', '#7e9af9', '#c084fc', '#f472b6'
        ]
        color_mapping = {sample: color_palette[i % len(color_palette)] for i, sample in enumerate(unique_samples)}
        color_mapping = {sample: color_mapping[sample.rsplit('_', 1)[0]] for sample in total_copy_number['Sample']}

        bar_chart = alt.Chart(total_copy_number).mark_bar(
            cornerRadius=8
        ).encode(
            x=alt.X('Sample', sort=x_sort, axis=alt.Axis(
                labelFontWeight='bold', 
                labelColor='#4B5563',
                titleFontWeight='bold', 
                titleColor='#374151',
                labelAngle=45
            )),
            y=alt.Y('Total Copy Number', axis=alt.Axis(
                labelFontWeight='bold',
                labelColor='#4B5563', 
                titleFontWeight='bold',
                titleColor='#374151'
            )),
            tooltip=['Sample', 'Total Copy Number'],
            color=alt.Color('Sample', scale=alt.Scale(
                domain=list(color_mapping.keys()), 
                range=list(color_mapping.values())
            ), legend=None)
        ).properties(
            width=400,
            height=400,
            title=alt.TitleParams(
                text='Total Copy Number per Sample',
                fontSize=28,
                fontWeight='bold',
                color='#1F2937'
            )
        )

        if region in st.session_state.pathogenic_TRs['region'].unique():
            pathogenic_threshold = st.session_state.pathogenic_TRs.loc[
                (st.session_state.pathogenic_TRs['region'] == region)
            ]['pathogenic_min'].values[0]

            threshold_line = alt.Chart(pd.DataFrame({'Total Copy Number': [pathogenic_threshold]})).mark_rule(
                color='#EF4444', strokeDash=[5, 5], size=2
            ).encode(y='Total Copy Number:Q')

            threshold_pointer = alt.Chart(pd.DataFrame({'Total Copy Number': [pathogenic_threshold]})).mark_text(
                text='🚨 Pathogenic Threshold', align='left', dx=5, dy=-10, fontSize=11, color='#DC2626', fontWeight='bold'
            ).encode(y='Total Copy Number:Q')
            
            bar_chart = alt.layer(bar_chart, threshold_line, threshold_pointer).configure_view(
                strokeWidth=0,
                fill='rgba(255,255,255,0.9)'
            )
        else:
            bar_chart = bar_chart.configure_view(
                strokeWidth=0,
                fill='rgba(255,255,255,0.9)'
            )

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
        

        return pd.DataFrame(data)

    def stack_plot(self, record, motif_names, sequences, span_list, motif_ids_list, sort_by="Value"):
        motif_colors = self.get_color_palette(len(record['motifs']))
        motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
        
        region = record['chr'] + ":" + str(record['pos'] - 1) + "-" + str(record['stop'] - 1)
        start = record['pos'] -1#+1
        stop = record['stop'] -1
        chrom = record['chr']
        updated_region = record['id']
        #st.write(updated_region)
        #updated_region = chrom + ":" + str(start) + "-" + str(stop)
        df = self.create_motif_dataframe(sequences, motif_colors, motif_ids_list, span_list, motif_names)
        if df.empty:
            return motif_colors, df
        df['Length'] = df['End'] - df['Start']
        
        default_height = 1200 
        chart_height = max(default_height, len(sequences) * 10)

        df['Sample'] = df['Sample'].apply(lambda x: x.replace("_pathogenic", ""))
        
        # Calculate total length per sample for sorting
        sample_lengths = df.groupby('Sample')['Length'].sum().sort_values(ascending=False)
        
        # Apply sorting based on the selected method
        if sort_by == "Value":
            # Sort by total length (descending)
            sorted_samples = sample_lengths.index.tolist()
        else:
            # Sort alphabetically
            sorted_samples = sorted(df['Sample'].unique())
        
        # Convert to categorical with the correct order
        df['Sample'] = pd.Categorical(df['Sample'], categories=sorted_samples, ordered=True)
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
        above_threshold_samples = None
        if updated_region in st.session_state.pathogenic_TRs['region'].unique():
            pathogenic_threshold = st.session_state.pathogenic_TRs.loc[
                (st.session_state.pathogenic_TRs['region'] == updated_region)
            ]['pathogenic_min'].values[0]
            if pathogenic_threshold is None:
                pathogenic_threshold = 0
            gene_name = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['gene'].values[0]
            inheritance = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['inheritance'].values[0]
            disease = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['disease'].values[0]

        # Display stats
        title = st.empty()


        title.markdown(f"""
            <div style="display: flex; gap: 20px; align-items: center; margin-bottom: 15px;">
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="font-size: 14px; color: #64748b; font-weight: 600;">Min:</div>
                    <div style="font-size: 16px; color: #1e293b; font-weight: 700; background: #f1f5f9; padding: 4px 12px; border-radius: 8px; border: 1px solid #e2e8f0;">{min_copy_number}</div>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="font-size: 14px; color: #64748b; font-weight: 600;">Max:</div>
                    <div style="font-size: 16px; color: #1e293b; font-weight: 700; background: #f1f5f9; padding: 4px 12px; border-radius: 8px; border: 1px solid #e2e8f0;">{max_copy_number}</div>
                </div>
            </div>
        """, unsafe_allow_html=True)
        # Prepare heatmap data
        heatmap_data = df[df['Motif'] != 'Interruption'].groupby(['Sample', 'Motif']).size().reset_index(name='Count')
        
        # Ensure heatmap data uses the same sample order as the main dataframe
        heatmap_data['Sample'] = pd.Categorical(heatmap_data['Sample'], categories=sorted_samples, ordered=True)
        heatmap_data = heatmap_data.sort_values('Sample')

        # Use the same sorting for both charts - based on the categorical order we already defined
        y_sort = sorted_samples  # This is the pre-defined order

        # Prepare stack plot data
        df["pathogenic"] = df["Total Copy Number"].apply(lambda x: "Pathogenic" if x >= pathogenic_threshold else "Not Pathogenic")
        pathogenic_threshold_length = pathogenic_threshold * len(motif_names[0])
        df['Sequence_length'] = df.groupby('Sample')['Length'].transform('sum')

        # Create heatmap with consistent sorting
        # Dynamically set the width of the heatmap based on the number of motifs
        motif_count = len(heatmap_data['Motif'].unique())
        # Base width per motif (adjust as needed, e.g., 50)
        width_per_motif = 40
        min_width = 40
        max_width = 480
        dynamic_width = max(min_width, min(width_per_motif * motif_count, max_width))

        heatmap = alt.Chart(heatmap_data).mark_rect(
            cornerRadius=4,
            stroke='white',
            strokeWidth=1
        ).encode(
            y=alt.Y('Sample:N', 
                title='', 
                sort=y_sort,  # Use the same sorting
                axis=alt.Axis(
                    labelFontWeight='bold', 
                    labelColor='#4B5563',
                    labelFontSize=20,
                    titleFontWeight='bold',
                    titleColor='#374151',
                    labelLimit=0,
                    ticks=False,
                    tickSize=20,
                    tickOffset=-2,
                    labelPadding=10,
                    labelOverlap=False,
                ),
                scale=alt.Scale(paddingInner=0, paddingOuter=0.1)),
            x=alt.X('Motif:N', 
                title='', 
                sort=alt.EncodingSortField(field='Count', op='sum', order='descending'),
                axis=alt.Axis(
                    labelFontWeight='bold',
                    labelFontSize=20,
                    labelColor='#4B5563',
                    titleFontWeight='bold',
                    titleColor='#374151',
                    labelAngle=90,
                    labelPadding=10,
                    labelLimit=100
                )),
            color=alt.Color('Count:Q', 
                        scale=alt.Scale(scheme='redpurple'), 
                        title='',
                        legend=alt.Legend(
                            title="",
                            orient="top",
                            titleFontSize=14,
                            titleColor='#1F2937',
                            gradientLength=120,
                            labelFontSize=12,
                            labelColor='#4B5563',
                            columns=4
                        )
            ),
            tooltip=['Sample', 'Motif', 'Count']
        ).properties(
            width=dynamic_width,
            height=chart_height,
            title=''
        )

        # Create stack plot with the same sorting
        stack_chart = alt.Chart(df).mark_bar(
            cornerRadius=0
        ).encode(
            y=alt.Y(
                'Sample:N',  # Use nominal type
                sort=y_sort,  # Use the same sorting
                title='',
                axis=alt.Axis(
                    labelOverlap=False, 
                    ticks=False, 
                    labelColor='#4B5563', 
                    labelFontSize=0, 
                    labelFontWeight='bold',
                    titleColor='#374151',
                    titleFontSize=28,
                    labelPadding = 1000,
                    labels = False,
                    labelLimit = 0,
                ),
                scale=alt.Scale(paddingInner=0, paddingOuter=0.6)
            ),
            x=alt.X('Length', 
                title='Sequence Length', 
                stack='zero', 
                axis=alt.Axis(
                    labelColor='#4B5563', 
                    labelFontSize=20, 
                    labelFontWeight='bold', 
                    titleColor='#374151',
                    titleFontSize=28,
                    labelOverlap=False,
                    labelAngle=90,  # Set labels to be 90 degrees
                    # Custom tick values: show every second number
                    tickMinStep=100,
                ), 
                scale=alt.Scale(domain=[0, df['Sequence_length'].max() + 10])),
            color=alt.Color(
                'Motif', 
                scale=alt.Scale(
                    domain=list(motif_names) + ['Interruption'], 
                    range=list(motif_colors.values()) + ['#EF4444']
                ),
                legend=alt.Legend(
                    title="",
                    orient="top",
                    gradientLength=120,
                    labelFontSize=20,
                    labelColor='#4B5563',
                    symbolStrokeWidth=3,
                    symbolSize=120,
                    symbolType="square",
                    columns=len(motif_names) + 1 if len(motif_names) + 1 <= 10 else 5,
                    # Increase padding between legend entries for more spacing
                    columnPadding=25
                )
            ),
            order=alt.Order('Order', sort='ascending'),
            tooltip=['Sample', 'Motif', 'Start', 'End', 'Sequence', 'pathogenic', 'Length', 'Sequence_length']
        ).properties(
            width=1300,
            height=chart_height,
            title='',
        )

        # Create pathogenic threshold elements
        threshold_elements = []
        
        # Always calculate above_threshold_samples, but only use it if threshold > 0
        above_threshold_samples = df.groupby('Sample')['Length'].sum().reset_index()
        if pathogenic_threshold > 0:
            above_threshold_samples = above_threshold_samples[above_threshold_samples['Length'] > pathogenic_threshold_length]
        else:
            above_threshold_samples = above_threshold_samples.iloc[0:0]  # Empty DataFrame with same structure

        if pathogenic_threshold > 0 and df.groupby('Sample')['Length'].sum().max() > pathogenic_threshold_length:
            # Add threshold line
            rule = alt.Chart(pd.DataFrame({'x': [pathogenic_threshold_length]})).mark_rule(
                color='#EF4444', strokeDash=[5, 5], strokeWidth=2
            ).encode(x='x:Q')
            threshold_elements.append(rule)
            
            if not above_threshold_samples.empty:
                # Create arrow data - position arrows slightly above each bar
                arrow_data = above_threshold_samples.copy()
                arrow_data['arrow_x'] = arrow_data['Length'] + (df['Sequence_length'].max() * 0.01)  # Small offset from bar end
                arrow_data['arrow_y_offset'] = 0.4  # Position in the middle of the bar vertically
                

                arrows = alt.Chart(arrow_data).mark_text(
                    text='➔',  # Right-arrow symbol
                    align='left',
                    baseline='middle',
                    fontSize=20,
                    color='#DC2626',
                    fontWeight='bold',
                    angle=235,   # Rotate the arrow to point right-up toward the bar
                    dy=0,
                    dx=-10        # Move the arrow a bit down
                ).encode(
                    x=alt.X('arrow_x:Q', title=None),
                    y=alt.Y('Sample:N', sort=y_sort),
                    tooltip=[alt.Tooltip('Sample:N'), alt.Tooltip('Length:Q', title='Total Length')]
                )
                
                threshold_elements.append(arrows)
                
                # Add warning text for pathogenic samples
                warning_text = alt.Chart(arrow_data).mark_text(
                    #text=' PATHOGENIC 🚨',
                    text=' PATHOGENIC',
                    align='left',
                    baseline='middle',
                    fontSize=12,
                    color='#DC2626',
                    fontWeight='bold',
                    dy=10,
                    dx=10,
                ).encode(
                    x=alt.X('arrow_x:Q', title=None),
                    y=alt.Y('Sample:N', sort=y_sort),
                    tooltip=[alt.Tooltip('Sample:N'), alt.Tooltip('Length:Q', title='Total Length')]
                )
                
                threshold_elements.append(warning_text)

        # Combine all elements
        if threshold_elements:
            # Add threshold elements to the stack chart
            all_stack_elements = [stack_chart] + threshold_elements
            final_stack_chart = alt.layer(*all_stack_elements)
        else:
            final_stack_chart = stack_chart

        # Create title with gene information if available
        chart_title = "Motif Distribution across Samples"
        subtitle_text = ""
        if gene_name or inheritance or disease or updated_region:
            gene_info_parts = []
            # In SVG exports, Unicode emoji/icons (like 🧬, 🏥, 🧭, etc.) may not render properly
            # Replace emoji with plain text or ASCII alternatives for SVG-export friendliness
            if gene_name:
                gene_info_parts.append(f"Gene: {gene_name}")
            if inheritance:
                gene_info_parts.append(f"Inheritance: {inheritance}")
            if disease:
                disease_display = disease #if len(disease) <= 50 else disease[:47] + "..."
                gene_info_parts.append(f"Disease: {disease_display}")
            if updated_region:
                gene_info_parts.append(f"Region: {updated_region}")
            subtitle_text = " • ".join(gene_info_parts)

        # Wrap the heatmap in a chart with a fixed width background if width is very small.
        min_heatmap_display_width = 120  # Adjust as desired for visual effect

        # If the calculated dynamic_width is exactly 40 (the minimum), pad the heatmap with a blank chart.
        if dynamic_width == 40:
            # Create a transparent dummy chart to the left of the heatmap to push it flush to the stack plot
            padding_width = min_heatmap_display_width - dynamic_width

            blank_left = alt.Chart(
                pd.DataFrame({"dummy": [0]})
            ).mark_rect(opacity=0).encode(
                x=alt.value(0),
                y=alt.value(0)
            ).properties(
                width=padding_width,
                height=chart_height
            )

            heatmap_display = alt.hconcat(
                blank_left,
                heatmap.properties(width=dynamic_width),
                spacing=0
            )
        else:
            heatmap_display = heatmap

        combined_chart = alt.hconcat(
            heatmap_display, 
            final_stack_chart,
            spacing=0  # Set gap between plots to 0 to keep them flush
        ).resolve_scale(
            y='shared'
        ).properties(
            title=alt.TitleParams(
                text=chart_title,
                fontSize=32,
                fontWeight='bold',
                anchor='middle',
                color='#1F2937',
                subtitle=subtitle_text,
                subtitleFontSize=20,
                subtitleColor='#4c1d95',
                subtitleFontWeight='bold',
                subtitlePadding=25,
            )
        ).configure_view(
            strokeWidth=0
        ).configure_scale(
            bandPaddingInner=0
        ).configure_axis(
            labelLimit=1000,
            tickCount=len(df['Sample'].unique()),
            labelOverlap=False,
            tickMinStep=max(200, int(df['Sequence_length'].max() * 0.1)),
        ).configure_axisY(
            labelFontSize=12,
            labelLimit=1000,
            tickCount=len(df['Sample'].unique()),
            labelOverlap=False,            
        )
        # Adjust for large number of samples
        if chart_height > default_height:
            combined_chart = combined_chart.configure_axisY(labelFontSize=8)
        
        # Display the combined chart
        st.altair_chart(combined_chart, use_container_width=True)
        
        # Show pathogenic info if applicable
        if pathogenic_threshold > 0 and above_threshold_samples is not None and not above_threshold_samples.empty:
            st.warning(f"🚨 **Pathogenic Alert**: {len(above_threshold_samples)} sample(s) exceed the pathogenic threshold of {pathogenic_threshold} copies")
        
        return motif_colors, df

    def count_motifs(self, motif_ids):
        motif_count = {}
        for idx, motif in enumerate(motif_ids):
            if motif in motif_count:
                motif_count[motif] += 1
            else:
                motif_count[motif] = 1
        return motif_count

    def visulize_TR_with_dynamic_sequence(self,record,hgsvc_records, container ,motif_colors,CN1_col,CN2_col, show_comparison):
        
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


    
 
        with container:
            st.markdown("""
                <style>
                    .stTabs [data-baseweb="tab-list"] {
                        gap: 6px;
                        background: #f8fafc;
                        padding: 4px 4px;
                        border-radius: 16px;
                        margin-bottom: 18px;
                        min-height: 22px;
                    }
                    .stTabs [data-baseweb="tab"] {
                        height: 28px;
                        min-height: 28px;
                        min-width: 38px;
                        white-space: pre-wrap;
                        background: white;
                        border-radius: 10px;
                        border: 2px solid #e2e8f0;
                        color: #64748b;
                        font-weight: 700;
                        font-size: 1.13rem;
                        padding: 0 13px;
                        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                        line-height: 1.33;
                    }
                    .stTabs [data-baseweb="tab"]:hover {
                        background: #f1f5f9;
                        border-color: #cbd5e1;
                        color: #475569;
                        transform: translateY(-1.8px) scale(1.04);
                        box-shadow: 0 2px 12px rgba(0,0,0,0.13);
                    }
                    .stTabs [aria-selected="true"] {
                        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                        color: white !important;
                        border-color: #667eea !important;
                        box-shadow: 0 5px 19px rgba(102, 126, 234, 0.23) !important;
                        transform: translateY(-2px) scale(1.045);
                        font-size: 1.23rem !important;
                    }
                    .stTabs [data-baseweb="tab-highlight"] {
                        background-color: transparent !important;
                    }
                </style>
            """, unsafe_allow_html=True)




        
            tab1, tab2, tab3 = st.tabs([
                "🧬 **Alleles**", 
                "🔄 **Alleles vs Ref**", 
                "🌍 **Alleles vs Pop**"
            ])
            with tab2:
                st.html('<div class="tab-content">')
                # Your tab1 content here
                motif_legend_html(record['motif_ids_ref'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("###  Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                display_dynamic_sequence_with_highlighted_motifs("Ref", record['ref_allele'], record['motif_ids_ref'], record['spans'][0], motif_colors, motif_names)
                alt_allele1 = record['alt_allele1']
                display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names, record['supported_reads_h1'])
                if record['alt_allele2'] != '':
                    alt_allele2 = record['alt_allele2'] 
                    display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names, record['supported_reads_h2'])
                st.html('</div>')
            with tab1:
                motif_legend_html(record['motif_ids_h1'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("###  Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                st.html('<div class="tab-content">')
                alt_allele1 = record['alt_allele1']
                display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names, record['supported_reads_h1'])

                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        with st.expander("Show motif bar plot for Allel1", expanded=False):
                            plot_motif_bar(motif_count_h1, motif_names, motif_colors)

                if record['alt_allele2'] != '':

                    alt_allele2 = record['alt_allele2'] 
                    display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names, record['supported_reads_h2'])
                    
                    plot_container_h2 = st.empty()
                    with plot_container_h2:
                        if show_comparison == False:
                            with st.expander("Show motif bar plot for Allel2", expanded=False):
                                plot_motif_bar(motif_count_h2, motif_names, motif_colors)
                st.html('</div>')
            with tab3:
                if hgsvc_records:
                    # Add genotype comparison for population data
                    st.markdown("### 🧬 Population Genotype Comparison")
                    population_genotypes = {"Current Sample": record['gt']}
                    for sample_name, sample_record in hgsvc_records.items():
                        if 'gt' in sample_record:
                            population_genotypes[sample_name] = sample_record['gt']
                    
                    if len(population_genotypes) > 1:
                        create_genotype_comparison_matrix(population_genotypes)
                        st.markdown("---")
                    
                    self.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
                else:
                    st.info("no population data found")

                st.html('</div>')
                
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


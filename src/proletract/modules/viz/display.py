"""
Display functions for visualization of tandem repeat regions and motifs.

This module provides functions for rendering genomic region displays, motif visualizations,
and interactive cohort/individual sample exploration interfaces.
"""
import streamlit as st
import pandas as pd
import re
from proletract.modules.viz import vis_helper as _vh
from proletract.modules.viz import parsers
from proletract.modules.viz import plots
from proletract.modules.viz import utils

# Re-export helper functions for convenience
motif_legend_html = _vh.motif_legend_html
plot_motif_bar = _vh.plot_motif_bar
display_motifs_as_bars = _vh.display_motifs_as_bars
display_genotype_card = _vh.display_genotype_card
display_dynamic_sequence_with_highlighted_motifs = _vh.display_dynamic_sequence_with_highlighted_motifs
create_genotype_comparison_matrix = _vh.create_genotype_comparison_matrix


def render_region_display(markdown_placeholder, region):
    """
    Render a genomic region display with external database links.
    
    Parses the region string to extract chromosome and position information,
    then generates clickable links to external genomic databases (UCSC, Ensembl,
    NCBI, gnomAD, DECIPHER, TRExplorer).
    
    Args:
        markdown_placeholder: Streamlit markdown placeholder for HTML injection.
        region: Genomic region string in format "chr:start-end" or "chr start end".
    """
    _vh.apply_global_typography()
    
    chrom = None
    pos_range = None
    chrom_no_chr = None
    start = None
    end = None
    
    try:
        if ':' in region:
            chrom, pos_range = region.split(':', 1)
            chrom_no_chr = chrom.replace('chr', '').replace('CHR', '')
            
            if '-' in pos_range:
                start, end = pos_range.split('-', 1)
                start = start.strip()
                end = end.strip()
            
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
    except Exception:
        pass
    
    urls = {}
    if chrom and pos_range:
        urls['UCSC'] = {
            'url': f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}:{pos_range}",
            'icon': '🌐',
            'color': 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)'
        }
        
        if chrom_no_chr:
            urls['Ensembl'] = {
                'url': f"https://www.ensembl.org/Homo_sapiens/Location/View?r={chrom_no_chr}:{pos_range}",
                'icon': '🧬',
                'color': 'linear-gradient(135deg, #FF6B6B 0%, #EE5A6F 100%)'
            }
            
            urls['NCBI'] = {
                'url': f"https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000001405.40&chr={chrom_no_chr}&from={start}&to={end}" if start and end else f"https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000001405.40&chr={chrom_no_chr}",
                'icon': '📊',
                'color': 'linear-gradient(135deg, #96CEB4 0%, #6C9A8B 100%)'
            }
            
            urls['DECIPHER'] = {
                'url': f"https://www.deciphergenomics.org/browser#q/grch37:{chrom_no_chr}:{pos_range}",
                'icon': '🔍',
                'color': 'linear-gradient(135deg, #A8E6CF 0%, #7FCDCD 100%)'
            }
        
        if chrom_no_chr and start and end:
            urls['gnomAD'] = {
                'url': f"https://gnomad.broadinstitute.org/region/{chrom_no_chr}-{start}-{end}",
                'icon': '📈',
                'color': 'linear-gradient(135deg, #FFE66D 0%, #FF6B6B 100%)'
            }
        
        if chrom and pos_range:
            urls['TRExplorer'] = {
                'url': f"https://trexplorer.broadinstitute.org/#sc=isPathogenic&sd=DESC&showRs=1&searchQuery={chrom}:{pos_range}&showColumns=0i1i2i3i4i7i21i17",
                'icon': '🔬',
                'color': 'linear-gradient(135deg, #10B981 0%, #059669 100%)'
            }
    
    html = f"""
        <style>
            .region-badge {{
                display: inline-block;
                background: #ECEAFB;
                padding: 4px 10px;
                border-radius: 36px;
                font-size: var(--pt-font-xs);
                font-weight: 600;
                margin-bottom: 14px;
                backdrop-filter: blur(8px);
                border: 1px solid #C9BEEF;
                letter-spacing: 0.04em;
            }}
            .external-link {{
                display: inline-flex;
                align-items: center;
                gap: 8px;
                color: white !important;
                text-decoration: none;
                padding: 12px 18px;
                border-radius: 14px;
                font-size: var(--pt-font-sm);
                font-weight: 600;
                transition: all 0.3s ease;
                box-shadow: 0 3px 10px rgba(0, 0, 0, 0.22);
            }}
            .external-link:hover {{
                transform: translateY(-2px) scale(1.02);
                box-shadow: 0 5px 14px rgba(0, 0, 0, 0.24);
            }}
            .external-link-icon {{
                font-size: var(--pt-font-md);
            }}
            .links-container {{
                display: flex;
                flex-wrap: wrap;
                gap: 14px;
                justify-content: center;
                margin-top: 16px;
            }}
        </style>
        
        <div style="display: flex; flex-direction: column; justify-content: center; align-items: center; min-height: 36px;">
            <div style="
                display: flex; 
                align-items: center; 
                background: #ECEAFB; 
                padding: 4px 10px; 
                border-radius: 12px; 
                box-shadow: 0 2px 6px rgba(118, 75, 162, 0.10); 
                font-size: var(--pt-font-lg); 
                font-weight: 600;
                border: 1px solid #C9BEEF;
                margin-bottom: 8px;
                ">
                <span style="width: 12px; height: 12px; background: #764ba2; border-radius: 50%; margin-right: 10px; box-shadow: 0 0 5px 1px #A184D6;"></span>
                Region: {region}
            </div>
            <div class="links-container" style="gap: 8px; margin-top: 8px;">
                {''.join([
                    f'<a href="{info["url"]}" target="_blank" class="external-link" style="background: {info["color"]}; padding: 6px 12px; border-radius: 8px; font-size: var(--pt-font-xs);">'
                    f'<span class="external-link-icon" style="font-size: var(--pt-font-sm);">{info["icon"]}</span> {name}</a>'
                    for name, info in urls.items()
                ])}
            </div>
        </div>
    """
    markdown_placeholder.html(html)


def display_motifs_with_bars(record, container, motif_colors, hgsvc_records):
    """
    Display tandem repeat motifs as bar visualizations.
    
    Renders motif counts for reference and alternate alleles using bar visualizations. 
    Args:
        record: VCF record dictionary containing motif and allele information.
        container: Streamlit container for rendering the visualization.
        motif_colors: Dictionary mapping motif indices to color codes.
        hgsvc_records: Dictionary of population reference records for comparison.
    """
    motif_names, motif_count_h1, motif_count_h2 = parsers.parse_motif_in_region(record, parsers.count_motifs)
    
    if motif_count_h1 is None and motif_count_h2 is None:
        st.info("No motifs found in the region")
        return
    
    # Use optimized CSS injection (cached, injected once)
    from .css_utils import inject_tab_styles_once
    inject_tab_styles_once()

    with container:
        
        tab1, tab2 = st.tabs([
            "🧬 **Alleles vs Ref**", 
            "🌍 **Alleles vs Pop**",
        ])
        # Tab content styling
        st.markdown(
            "<style>.tab-content {font-size: var(--pt-font-md);}</style>",
            unsafe_allow_html=True,
        )

        with tab1:
            motif_legend_html(record['motif_ids_ref'], motif_colors, motif_names)
            display_genotype_card(record['gt'], "Current Sample", show_details=True)
            
            display_motifs_as_bars("Ref", motif_colors, record['motif_ids_ref'], record['spans'][0], record['ref_allele'], motif_names, None)
            display_motifs_as_bars("Allel1", motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names, record['supported_reads_h1'])
            if record['alt_allele2'] != '':
                display_motifs_as_bars("Allel2", motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names, record['supported_reads_h2'])

        with tab2:
            # Population comparison tab content
            if hgsvc_records:
                current_sample_gt = record['gt']
                # Split the genotype into haplotype genotypes
                gt_parts = current_sample_gt.split('/') if '/' in str(current_sample_gt) else [current_sample_gt, current_sample_gt]
                gt_h1 = gt_parts[0] if len(gt_parts) > 0 else "0"
                gt_h2 = gt_parts[1] if len(gt_parts) > 1 else gt_parts[0]
                
                # Compute the diploid genotype assembly
                diploid_gt = parsers.compute_diploid_genotype_assembly(
                    gt_h1, gt_h2, 
                    record.get('motif_ids_h1', []), 
                    record.get('motif_ids_h2', [])
                )
                
                # Create a dictionary of population genotypes
                population_genotypes = {"Current Sample": diploid_gt}
                
                # Create a dictionary of sample groups
                sample_groups = {}
                # Iterate over the hgsvc records
                for sample_name, sample_record in hgsvc_records.items():
                    if sample_name.endswith('_h1'):
                        base_name = sample_name[:-3]
                        # Create a dictionary of sample groups
                        if base_name not in sample_groups:
                            sample_groups[base_name] = {}
                        sample_groups[base_name]['h1'] = sample_record
                    elif sample_name.endswith('_h2'):
                        base_name = sample_name[:-3]
                        if base_name not in sample_groups:
                            sample_groups[base_name] = {}
                        sample_groups[base_name]['h2'] = sample_record
                    else:
                        if 'gt' in sample_record:
                            population_genotypes[sample_name] = sample_record['gt']
                
                # Iterate over the sample groups
                for base_name, haplotypes in sample_groups.items():
                    # If both haplotypes are present
                    if 'h1' in haplotypes and 'h2' in haplotypes:
                        # Get the h1 record
                        h1_record = haplotypes['h1']
                        # Get the h2 record
                        h2_record = haplotypes['h2']
                        
                        # Get the genotype of the h1 record
                        gt_h1_pop = h1_record.get('gt', '0')
                        # Get the genotype of the h2 record
                        gt_h2_pop = h2_record.get('gt', '0')
                        # Get the motif ids of the h1 record
                        ids_h1_pop = h1_record.get('motif_ids_h', [])
                        # Get the motif ids of the h2 record
                        ids_h2_pop = h2_record.get('motif_ids_h', [])
                        
                        # Compute the diploid genotype assembly
                        diploid_gt_pop = parsers.compute_diploid_genotype_assembly(gt_h1_pop, gt_h2_pop, ids_h1_pop, ids_h2_pop)
                        # Add the diploid genotype to the population genotypes
                        population_genotypes[base_name] = diploid_gt_pop
                    elif 'h1' in haplotypes:
                        # If only the h1 haplotype is present
                        population_genotypes[base_name] = haplotypes['h1'].get('gt', '0')
                        # If only the h2 haplotype is present
                    elif 'h2' in haplotypes:
                        # If only the h2 haplotype is present
                        population_genotypes[base_name] = haplotypes['h2'].get('gt', '0')
                
                # If there are more than one population genotypes
                if len(population_genotypes) > 1:
                    # Create a genotype comparison matrix
                    create_genotype_comparison_matrix(population_genotypes)
                    # Add a separator
                    st.markdown("---")
                
                # Plot the HGSVC vs allele
                plots.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
            else:
                st.info("no population data found")


def display_motif_legend(motifs, motif_colors, right_column):
    """
    Display a scrollable legend showing motif colors and names.
    
    Args:
        motifs: List or tuple of motif sequences.
        motif_colors: Dictionary mapping motif indices to color codes.
        right_column: Streamlit column container for the legend.
    """
    st.markdown("### Motif Legend")
    st.markdown('<div style="max-height:400px; overflow-y:scroll;">', unsafe_allow_html=True)
    
    if isinstance(motifs, tuple):
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


def visulize_TR_with_dynamic_sequence(record, hgsvc_records, container, motif_colors, CN1_col, CN2_col, show_comparison):
    """
    Visualize tandem repeat regions with dynamically highlighted motif sequences.
    
    Displays reference and alternate alleles with color-coded motifs in a tabbed
    interface, with options for population comparison.
    
    Args:
        record: VCF record dictionary containing sequence and motif information.
        hgsvc_records: Dictionary of population reference records.
        container: Streamlit container for rendering.
        motif_colors: Dictionary mapping motif indices to color codes.
        CN1_col: Streamlit column for copy number display (unused).
        CN2_col: Streamlit column for copy number display (unused).
        show_comparison: Boolean flag for comparison mode (unused).
    """
    if record['motif_ids_h1'] == ['.'] and record['motif_ids_h2'] == ['.']:
        st.info("No motifs found in the region")
        return
    
    motif_names = record['motifs']
    reference_copy_number = record['ref_CN']
    motif_count_ref = parsers.count_motifs(record['motif_ids_ref'])
    found_motifs_ref = list(motif_count_ref.keys())
    found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]

    motif_count_h1 = parsers.count_motifs(record['motif_ids_h1'])
    found_motifs_h1 = list(motif_count_h1.keys())
    found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
    motif_count_h2 = parsers.count_motifs(record['motif_ids_h2'])
    found_motifs_h2 = list(motif_count_h2.keys())
    found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
    motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
    total_copy_number_h1 = str(record['spans'][1]).count('-')

    with container:
        from .css_utils import inject_tab_styles_once
        inject_tab_styles_once()

        tab1, tab2 = st.tabs([
            "🧬 **Alleles vs Ref**", 
            "🌍 **Alleles vs Pop**", 
        ])
        
        with tab1:
            st.html('<div class="tab-content">')
            motif_legend_html(record['motif_ids_ref'], motif_colors, motif_names)
            display_genotype_card(record['gt'], "Current Sample", show_details=True)
            
            display_dynamic_sequence_with_highlighted_motifs("Ref", record['ref_allele'], record['motif_ids_ref'], record['spans'][0], motif_colors, motif_names)
            alt_allele1 = record['alt_allele1']
            display_dynamic_sequence_with_highlighted_motifs("Allel1", alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names, record['supported_reads_h1'])
            if record['alt_allele2'] != '':
                alt_allele2 = record['alt_allele2']
                display_dynamic_sequence_with_highlighted_motifs("Allel2", alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names, record['supported_reads_h2'])
            st.html('</div>')

        with tab2:
            if hgsvc_records:
                current_sample_gt = record['gt']
                gt_parts = current_sample_gt.split('/') if '/' in str(current_sample_gt) else [current_sample_gt, current_sample_gt]
                gt_h1 = gt_parts[0] if len(gt_parts) > 0 else "0"
                gt_h2 = gt_parts[1] if len(gt_parts) > 1 else gt_parts[0]
                
                diploid_gt = parsers.compute_diploid_genotype_assembly(
                    gt_h1, gt_h2, 
                    record.get('motif_ids_h1', []), 
                    record.get('motif_ids_h2', [])
                )
                
                population_genotypes = {"Current Sample": diploid_gt}
                
                sample_groups = {}
                for sample_name, sample_record in hgsvc_records.items():
                    if sample_name.endswith('_h1'):
                        base_name = sample_name[:-3]
                        if base_name not in sample_groups:
                            sample_groups[base_name] = {}
                        sample_groups[base_name]['h1'] = sample_record
                    elif sample_name.endswith('_h2'):
                        base_name = sample_name[:-3]
                        if base_name not in sample_groups:
                            sample_groups[base_name] = {}
                        sample_groups[base_name]['h2'] = sample_record
                    else:
                        if 'gt' in sample_record:
                            population_genotypes[sample_name] = sample_record['gt']
                
                for base_name, haplotypes in sample_groups.items():
                    if 'h1' in haplotypes and 'h2' in haplotypes:
                        h1_record = haplotypes['h1']
                        h2_record = haplotypes['h2']
                        
                        gt_h1_pop = h1_record.get('gt', '0')
                        gt_h2_pop = h2_record.get('gt', '0')
                        ids_h1_pop = h1_record.get('motif_ids_h', [])
                        ids_h2_pop = h2_record.get('motif_ids_h', [])
                        
                        diploid_gt_pop = parsers.compute_diploid_genotype_assembly(gt_h1_pop, gt_h2_pop, ids_h1_pop, ids_h2_pop)
                        population_genotypes[base_name] = diploid_gt_pop
                    elif 'h1' in haplotypes:
                        population_genotypes[base_name] = haplotypes['h1'].get('gt', '0')
                    elif 'h2' in haplotypes:
                        population_genotypes[base_name] = haplotypes['h2'].get('gt', '0')
                
                if len(population_genotypes) > 1:
                    create_genotype_comparison_matrix(population_genotypes)
                    st.markdown("---")
                
                plots.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
            else:
                st.info("no population data found")

            st.html('</div>')


def visulize_cohort():
    """
    Main visualization function for cohort mode.
    
    Provides an interactive interface for exploring tandem repeat regions across
    multiple samples, with search functionality and navigation controls.
    """
    if 'cohorts_records_map' not in st.session_state:
        return
    
    region_options = list(st.session_state.cohorts_records_map.values())
    
    regions_idx = st.session_state.get('regions_idx', 0)
    if regions_idx >= len(region_options):
        regions_idx = 0
        st.session_state.regions_idx = 0
    
    default_region = region_options[regions_idx] if region_options else ""
    
    if 'cached_region_options_cohort' not in st.session_state:
        st.session_state.cached_region_options_cohort = region_options
    
    if 'region_selected_cohort' not in st.session_state:
        st.session_state.region_selected_cohort = ""
    
    search_query = st.sidebar.text_input(
        "🔍 Search region:", 
        value=st.session_state.region_selected_cohort,
        key="region_search_cohort",
        help="Type to search and filter results",
        placeholder="Type to search..."
    )
    
    st.markdown("""
        <style>
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
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] input {
                color: black !important;
                font-weight: 700;
            }
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[aria-selected="true"] span {
                color: black !important;
            }
        </style>
    """, unsafe_allow_html=True)

    if search_query:
        search_lower = search_query.lower()
        filtered = [r for r in st.session_state.cached_region_options_cohort if search_lower in r.lower()]
        filtered_regions = filtered[:10]
        total_matches = len(filtered)
        
        if filtered_regions:
            region = st.sidebar.selectbox(
                " ",
                filtered_regions, 
                index=0,
                key="region_suggest",
                help="Select from suggestions",
                label_visibility="collapsed"
            )
            st.session_state.region_selected_cohort = region
        else:
            region = search_query
            
        if total_matches > 10:
            st.sidebar.markdown(f"<span style='font-size:11px; color:orange;'>Showing 10 of {total_matches:,} matches</span>", unsafe_allow_html=True)
        else:
            st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>{total_matches:,} matches</span>", unsafe_allow_html=True)
    else:
        region = default_region
        total_matches = len(st.session_state.cached_region_options_cohort)
        st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>Search to find from {total_matches:,} regions</span>", unsafe_allow_html=True)

    if 'regions_idx' not in st.session_state:
        st.session_state.regions_idx = 0
    
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
            st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.cohorts_records_map) - 1)
    
    st.markdown("""
        <style>
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
        region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
    
    mode_placeholder = st.empty()
    
    middel, _ = st.columns([1, 0.1], gap="small")
    render_region_display(middel, region)
    st.markdown("""
        <style>
            :root {
                --region-color-light: black;
                --region-color-dark: white;
            }
            .region-container {
                color: var(--region-color-light);
            }
            @media (prefers-color-scheme: dark) {
                .region-container {
                    color: var(--region-color);
                }
            }
        </style>
    """, unsafe_allow_html=True)

    st.session_state.cohort_results = utils.get_results_cohort(
        region, st.session_state.cohort_files, st.session_state.cohort_file_paths, 
        st.session_state.cohort_mode,
        parsers.parse_record,
        parsers.parse_record_assembly
    )
    
    if 'cohort_results' in st.session_state:
        region = st.session_state.regions_idx
        if st.session_state.cohort_results == {}:
            st.warning(f"No motifs found in the region: {region}")
            st.stop()
        else:
            plots.plot_Cohort_results(st.session_state.cohort_results, st.session_state.cohort_mode)
    else:
        st.stop()


def visulize_region():
    """
    Main visualization function for individual sample mode.
    
    Provides an interactive interface for exploring tandem repeat regions in a
    single sample, with search functionality, navigation controls, and multiple
    display options (sequence view or bar view).
    """
    if 'regions_idx' not in st.session_state:
        st.session_state.regions_idx = 0
    
    region_options = list(st.session_state.records_map.values())
    
    regions_idx = st.session_state.get('regions_idx', 0)
    if regions_idx >= len(region_options):
        regions_idx = 0
        st.session_state.regions_idx = 0
    
    default_region = region_options[regions_idx] if region_options else ""
    st.sidebar.markdown("### Select Region to Visualize")
    
    if 'cached_region_options' not in st.session_state:
        st.session_state.cached_region_options = region_options
    
    if 'region_selected_ind' not in st.session_state:
        st.session_state.region_selected_ind = ""
    
    search_query = st.sidebar.text_input(
        "🔍 Search region:", 
        value=st.session_state.region_selected_ind,
        key="region_search",
        help="Type to search and filter results",
        placeholder="Type to search..."
    )
    
    if search_query:
        search_lower = search_query.lower()
        filtered = [r for r in st.session_state.cached_region_options if search_lower in r.lower()]
        filtered_regions = filtered[:10]
        total_matches = len(filtered)
        
        if filtered_regions:
            region = st.sidebar.selectbox(
                " ",
                filtered_regions, 
                index=0,
                key="region_suggest",
                help="Select from filtered results",
                label_visibility="collapsed"
            )
            st.session_state.region_selected_ind = region
        else:
            region = search_query
            
        if total_matches > 10:
            st.sidebar.markdown(f"<span style='font-size:11px; color:orange;'>Showing 10 of {total_matches:,} matches</span>", unsafe_allow_html=True)
        else:
            st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>{total_matches:,} matches</span>", unsafe_allow_html=True)
    else:
        region = default_region
        total_matches = len(st.session_state.cached_region_options)
        st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>Search to find from {total_matches:,} regions</span>", unsafe_allow_html=True)
    
    display_option = st.sidebar.radio("Select Display Type", 
                                      ("Sequence with Highlighted Motifs", "Bars"))

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
    
    st.markdown("""
        <style>
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
    record = parsers.parse_record(st.session_state.vcf_file_path, record_key)
    
    if record is None:
        st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
        st.stop()
    
    hgsvc_records = utils.get_results_hgsvc_pop(record_key, st.session_state.files, st.session_state.file_paths, parsers.parse_record_assembly)
    
    if len(record["motif_ids_h1"]) == 0 and len(record["motif_ids_h2"]) == 0:
        st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
        st.stop()

    render_region_display(middel, st.session_state.records_map[st.session_state.regions_idx])

    st.markdown("""
        <style>
            :root {
                --region-color-light: black;
                --region-color-dark: white;
            }
            .region-container {
                color: var(--region-color-light);
            }
            @media (prefers-color-scheme: dark) {
                .region-container {
                    color: var(--region-color);
                }
            }
        </style>
    """, unsafe_allow_html=True)

    container = st.container()
    motif_colors = utils.get_color_palette(len(record['motifs']))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}

    if display_option == "Sequence with Highlighted Motifs":
        visulize_TR_with_dynamic_sequence(record, hgsvc_records, container, motif_colors, CN1_col, CN2_col, st.session_state.get('show_comparison', False))
    elif display_option == "Bars":
        display_motifs_with_bars(record, container, motif_colors, hgsvc_records)

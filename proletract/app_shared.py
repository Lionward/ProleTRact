import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path
import pandas as pd
import pysam
import plotly.express as px
import plotly.graph_objects as go
from proletract.modules.io.data_handling import VCFHandler, CohortHandler
from proletract.modules.viz.visualization import Visualization


def get_statistics_color_palette():
    """Get color palette for statistics plots matching app theme"""
    return [
        '#667eea',  # Primary blue
        '#764ba2',  # Primary purple
        '#4fd1c7',  # Teal
        '#68d391',  # Green
        '#f6ad55',  # Orange
        '#fc8181',  # Coral
        '#f093fb',  # Pink
        '#f6e05e',  # Yellow
        '#5a67d8',  # Darker blue
        '#6b46c1',  # Darker purple
        '#38a169',  # Darker green
        '#dd6b20',  # Darker orange
        '#e53e3e',  # Darker red
        '#d53f8c',  # Darker pink
        '#d69e2e',  # Darker yellow
        '#319795',  # Darker teal
    ]


def ensure_state_initialized():
    if 'vcf_handler' not in st.session_state:
        st.session_state.vcf_handler = VCFHandler()
    if 'cohort_handler' not in st.session_state:
        st.session_state.cohort_handler = CohortHandler()
    st.session_state.visualization = Visualization()
    get_pathogenic_TRs()


def get_pathogenic_TRs():
    root = Path(__file__).resolve().parent.parent
    data_path = root / "proletract" / "data" / "pathogenic_TRs.bed"
    pathogenic_trs = pd.read_csv(str(data_path), sep="\t", header=None)
    pathogenic_trs.columns = ["chrom", 'start', 'end', 'motif', 'pathogenic_min', 'inheritance', 'disease', 'gene']
    st.session_state.pathogenic_TRs = pathogenic_trs
    pathogenic_trs['region'] = pathogenic_trs['chrom'] + ":" + pathogenic_trs['start'].astype(str) + "-" + pathogenic_trs['end'].astype(str)


def configure_page():
    st.set_page_config(layout="wide")
    inject_global_styles()
    render_sidebar_branding()


def inject_global_styles():
    st.markdown("""
        <style>
            /* Slightly increase global font sizes */
            html, body { font-size: 24px; }
            [data-testid="stAppViewContainer"] h1 { font-size: 2.15rem !important; }
            [data-testid="stAppViewContainer"] h2 { font-size: 1.7rem !important; }
            [data-testid="stAppViewContainer"] h3 { font-size: 1.35rem !important; }
            [data-testid="stAppViewContainer"] p, 
            [data-testid="stAppViewContainer"] li, 
            [data-testid="stAppViewContainer"] code, 
            [data-testid="stAppViewContainer"] .stMarkdown, 
            [data-testid="stAppViewContainer"] label { 
                font-size: 1.05rem !important; 
            }
            /* Ensure visualization-specific components scale with global font */
            [data-testid="stAppViewContainer"] .sequence-dashboard, 
            [data-testid="stAppViewContainer"] .sequence-dashboard * {
                font-size: 1rem !important;
            }
            [data-testid="stAppViewContainer"] .sequence-header {
                font-size: 1.05rem !important;
            }
            [data-testid="stAppViewContainer"] .sequence-length {
                font-size: 0.98rem !important;
            }
            [data-testid="stAppViewContainer"] .motif-legend-container,
            [data-testid="stAppViewContainer"] .motif-legend-container * {
                font-size: 1rem !important;
            }
            [data-testid="stAppViewContainer"] .legend-stats .stat-item span,
            [data-testid="stAppViewContainer"] .legend-stats .stat-item .stat-value {
                font-size: 1rem !important;
            }
            /* Dataframe/table text */
            [data-testid="stAppViewContainer"] [data-testid="stDataFrame"] * { 
                font-size: 0.98rem !important; 
            }
            /* Sidebar base font size */
            [data-testid="stSidebar"] { font-size: 1.06rem !important; }
            [data-testid="stSidebar"] h1 { font-size: 1.9rem !important; }
            [data-testid="stSidebar"] h2 { font-size: 1.5rem !important; }
            [data-testid="stSidebar"] h3 { font-size: 1.25rem !important; }
            div[data-testid="stAppViewContainer"] .stButton > button {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                border: none;
                padding: 8px 14px;
                border-radius: 10px;
                font-weight: 700 !important;
                font-size: 16px !important;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
                margin: 4px;
                min-width: 130px;
                min-height: 44px;
                letter-spacing: 0.65px;
                line-height: 1.15 !important;
            }
            div[data-testid="stAppViewContainer"] .stButton > button * {
                font-size: 18px !important;
                font-weight: 700 !important;
                line-height: 1.15 !important;
            }
            div[data-testid="stAppViewContainer"] .stButton > button:hover {
                transform: translateY(-2px) scale(1.03);
                box-shadow: 0 8px 18px rgba(102, 126, 234, 0.28);
                background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
                font-size: 17px !important;
            }
            div[data-testid="stAppViewContainer"] .stButton > button:hover * {
                font-size: 17px !important;
            }
            div[data-testid="stAppViewContainer"] .stButton > button:active,
            div[data-testid="stAppViewContainer"] .stButton > button:focus {
                transform: scale(0.97);
                box-shadow: 0 2px 8px rgba(102, 126, 234, 0.18);
                font-size: 16px !important;
            }
            div[data-testid="stAppViewContainer"] .stButton > button:active *,
            div[data-testid="stAppViewContainer"] .stButton > button:focus * {
                font-size: 16px !important;
            }
            div[data-testid="stAppViewContainer"] .stButton button > div,
            div[data-testid="stAppViewContainer"] .stButton button > span,
            div[data-testid="stAppViewContainer"] .stButton button > div > span {
                font-size: 16px !important;
                font-weight: 700 !important;
                line-height: 1.15 !important;
            }
            section[data-testid="stSidebar"] .stButton > button {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                border: none;
                padding: 12px 18px;
                border-radius: 10px;
                font-weight: 600 !important;
                font-size: 18px !important;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                box-shadow: 0 2px 6px rgba(102, 126, 234, 0.12);
                margin: 3px 0;
                min-width: 140px;
                min-height: 48px;
                letter-spacing: 0.3px;
                line-height: 1.1 !important;
            }
            section[data-testid="stSidebar"] .stButton > button * {
                font-size: 18px !important;
                font-weight: 600 !important;
                line-height: 1.1 !important;
            }
            section[data-testid="stSidebar"] .stButton > button:hover {
                transform: translateY(-1.5px) scale(1.01);
                box-shadow: 0 3px 8px rgba(102, 126, 234, 0.15);
                background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
                font-size: 19px !important;
            }
            section[data-testid="stSidebar"] .stButton > button:hover * {
                font-size: 19px !important;
            }
            section[data-testid="stSidebar"] .stButton > button:active,
            section[data-testid="stSidebar"] .stButton > button:focus {
                transform: scale(0.97);
                box-shadow: 0 1px 3px rgba(102, 126, 234, 0.10);
                font-size: 17px !important;
            }
            section[data-testid="stSidebar"] .stButton > button:active *,
            section[data-testid="stSidebar"] .stButton > button:focus * {
                font-size: 17px !important;
            }
            [data-testid="stSidebar"] { background: linear-gradient(180deg, #667eea 0%, #764ba2 100%); }
            [data-testid="stSidebar"] * { color: white; }
            /* Ensure typed text in sidebar inputs is dark/legible */
            [data-testid="stSidebar"] input,
            [data-testid="stSidebar"] textarea,
            [data-testid="stSidebar"] select,
            [data-testid="stSidebar"] .stTextInput input,
            [data-testid="stSidebar"] .stTextArea textarea,
            [data-testid="stSidebar"] .stSelectbox div[role="combobox"],
            [data-testid="stSidebar"] .stMultiSelect div[role="combobox"],
            [data-testid="stSidebar"] .stNumberInput input {
                color: #1f2937 !important; /* dark gray */
                background: rgba(255, 255, 255, 0.96) !important;
            }
            [data-testid="stSidebar"] input::placeholder,
            [data-testid="stSidebar"] textarea::placeholder {
                color: #6b7280 !important; /* slate-500 */
            }
            /* Labels inside white input containers should be dark */
            [data-testid="stSidebar"] label,
            [data-testid="stSidebar"] legend,
            [data-testid="stSidebar"] .stMarkdown p,
            [data-testid="stSidebar"] .stCheckbox span,
            [data-testid="stSidebar"] .stRadio span {
                color: #1f2937 !important;
            }
            /* File uploader text colors */
            [data-testid="stSidebar"] [data-testid="stFileUploader"] *,
            [data-testid="stSidebar"] [data-testid="stFileUploader"] label {
                color: #1f2937 !important;
            }
            /* Hide Streamlit's default multipage navigation to control ordering */
            [data-testid="stSidebarNav"] { display: none !important; }
            section[data-testid="stSidebar"] nav { display: none !important; }
        </style>
    """, unsafe_allow_html=True)

    components.html("""
        <style>
        .floating-nav-btn {
            position: fixed !important;
            bottom: 24px;
            z-index: 9999 !important;
            min-width: 92px;
            min-height: 72px;
            background: linear-gradient(108deg, #764ba2 6%, #667eea 94%);
            color: #fff !important;
            border: none !important;
            border-radius: 22px !important;
            box-shadow: 0 8px 36px 0 rgba(138,105,227,0.25);
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            font-weight: 900;
            font-size: 19px;
            letter-spacing: 1.1px;
            cursor: pointer;
            opacity: 0.992;
            padding: 0 24px 12px 24px;
            gap: 3px;
            transition: box-shadow 0.2s, background 0.2s, transform 0.14s;
        }
        .floating-nav-btn:hover {
            background: linear-gradient(111deg, #667eea 8%, #764ba2 92%);
            transform: translateY(-8px) scale(1.07) rotate(-2.5deg);
        }
        </style>
        <script>
        (function() {
            const buttons = Array.from(document.querySelectorAll('button'));
            let prevBtn, nextBtn;
            buttons.forEach(btn => {
                const tx = btn.textContent || btn.innerText || '';
                if (tx.includes('Previous') && !btn.classList.contains('floating-nav-btn')) prevBtn = btn;
                if (tx.includes('Next') && !btn.classList.contains('floating-nav-btn')) nextBtn = btn;
            });
            [prevBtn, nextBtn].forEach((btn, i) => {
                if (btn) {
                    btn.classList.add('floating-nav-btn');
                    btn.style.left = i === 0 ? 'calc(50% - 116px)' : 'calc(50% + 44px)';
                    btn.style.right = '';
                    btn.style.position = 'fixed';
                    btn.style.bottom = '26px';
                }
            });
        })();
        </script>
    """, height=0)


@st.cache_data
def compute_vcf_statistics(vcf_file_path):
    """Compute statistics from VCF file (cached)"""
    if not vcf_file_path or not Path(vcf_file_path).exists():
        return None
    
    try:
        vcf = pysam.VariantFile(vcf_file_path)
        
        # Collect statistics
        total_regions = 0
        motif_size_counts = {}  # max_motif_size -> count
        motif_lengths = []  # all motif lengths for distribution
        regions_by_chromosome = {}
        genotype_counts = {}  # genotype -> count
        
        for rec in vcf.fetch():
            total_regions += 1
            
            # Get motifs from INFO field
            motifs = rec.info.get('MOTIFS', [])
            if isinstance(motifs, tuple):
                motifs = list(motifs)
            elif not isinstance(motifs, list):
                motifs = [motifs] if motifs else []
            
            # Calculate max motif size for this region
            max_motif_size = 0
            if motifs:
                motif_lengths_in_region = [len(str(m)) for m in motifs if m]
                if motif_lengths_in_region:
                    max_motif_size = max(motif_lengths_in_region)
                    motif_lengths.extend(motif_lengths_in_region)
            
            # Categorize by max motif size
            if max_motif_size == 0:
                category = "Unknown"
            elif max_motif_size == 1:
                category = "1"
            elif max_motif_size == 2:
                category = "2"
            elif max_motif_size == 3:
                category = "3"
            elif max_motif_size == 4:
                category = "4"
            elif max_motif_size == 5:
                category = "5"
            elif max_motif_size == 6:
                category = "6"
            elif max_motif_size == 7:
                category = "7"
            elif max_motif_size == 8:
                category = "8"
            elif max_motif_size == 9:
                category = "9"
            elif max_motif_size == 10:
                category = "10"
            else:
                category = ">10"
            
            motif_size_counts[category] = motif_size_counts.get(category, 0) + 1
            
            # Count by chromosome
            chrom = rec.chrom
            regions_by_chromosome[chrom] = regions_by_chromosome.get(chrom, 0) + 1
            
            # Extract genotype information
            if len(rec.samples) > 0 and 'GT' in rec.samples[0]:
                gt = rec.samples[0]['GT']
                if gt is not None:
                    # Handle different GT formats
                    if isinstance(gt, tuple):
                        gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                    elif isinstance(gt, (list,)):
                        gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                    else:
                        gt_str = str(gt)
                    genotype_counts[gt_str] = genotype_counts.get(gt_str, 0) + 1
        
        return {
            'total_regions': total_regions,
            'motif_size_counts': motif_size_counts,
            'motif_lengths': motif_lengths,
            'regions_by_chromosome': regions_by_chromosome,
            'genotype_counts': genotype_counts
        }
    except Exception as e:
        return {'error': str(e)}


def display_vcf_statistics(vcf_file_path, mode='individual'):
    """Display statistics and plots about TR regions from loaded VCF file"""
    stats = compute_vcf_statistics(vcf_file_path)
    
    if stats is None:
        st.info("Please load a VCF file to view statistics.")
        return
    
    if 'error' in stats:
        st.error(f"Error computing statistics: {stats['error']}")
        return
    
    total_regions = stats['total_regions']
    motif_size_counts = stats['motif_size_counts']
    motif_lengths = stats['motif_lengths']
    regions_by_chromosome = stats['regions_by_chromosome']
    genotype_counts = stats.get('genotype_counts', {})
    
    if total_regions == 0:
        st.warning("No regions found in VCF file.")
        return
    
    # Display statistics
    st.markdown(
        """
        <div style="border-radius: 0.5rem; border: 2px solid #667eea; background: #f7faff; padding: 1.2rem 1.5rem; margin-bottom:0.5rem;">
            <span style="font-size: 1.8rem; font-weight:700; color: #243c5a;">📊 General Statistics</span>
        </div>
        """,
        unsafe_allow_html=True
    )
    # Summary metrics - Row 1 with larger fonts
    st.markdown("""
        <style>
        [data-testid="stMetricValue"] {
            font-size: 2rem !important;
        }
        [data-testid="stMetricLabel"] {
            font-size: 1.2rem !important;
        }
        </style>
    """, unsafe_allow_html=True)
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Regions", f"{total_regions:,}")
    with col2:
        st.metric("Chromosomes", len(regions_by_chromosome))
    with col3:
        if motif_lengths:
            avg_motif_size = sum(motif_lengths) / len(motif_lengths)
            st.metric("Avg Motif Size", f"{avg_motif_size:.1f} bp")
        else:
            st.metric("Avg Motif Size", "N/A")
    with col4:
        if motif_lengths:
            max_overall = max(motif_lengths)
            st.metric("Max Motif Size", f"{max_overall} bp")
        else:
            st.metric("Max Motif Size", "N/A")
    
    # Genotype statistics - Row 2
    if genotype_counts:
        with st.container():
            st.markdown(
                """
                <div style="border-radius: 0.5rem; border: 2px solid #764ba2; background: #faf6ff; padding: 1.2rem 1.5rem; margin-bottom:0.5rem;">
                    <span style="font-size: 1.75rem; font-weight: 700; color: #4B247B;">🧬 Genotype Distribution</span>
                </div>
                """,
                unsafe_allow_html=True
            )
        # Sort genotypes logically (0/0, 0/1, 1/1, etc.)
        sorted_genotypes = sorted(genotype_counts.items(), 
                                 key=lambda x: (len(x[0].split('/')), x[0]))
        
        # Create columns for genotype metrics (up to 6 most common)
        num_genotype_cols = min(len(sorted_genotypes), 6)
        if num_genotype_cols > 0:
            genotype_cols = st.columns(num_genotype_cols)
            for idx, (gt, count) in enumerate(sorted_genotypes[:num_genotype_cols]):
                with genotype_cols[idx]:
                    percentage = (count / total_regions) * 100
                    st.metric(f"Genotype {gt}", f"{count:,}", f"{percentage:.1f}%")
            
            # If there are more genotypes, show them in a dataframe
            if len(sorted_genotypes) > 6:
                remaining_genotypes = sorted_genotypes[6:]
                st.markdown("""
                    <div style="font-size: 1.25rem; font-weight: 600; margin-top: 1rem; margin-bottom: 0.5rem;">
                        <strong>Other genotypes:</strong>
                    </div>
                """, unsafe_allow_html=True)
                df_other_gt = pd.DataFrame({
                    'Genotype': [gt for gt, _ in remaining_genotypes],
                    'Count': [count for _, count in remaining_genotypes],
                    'Percentage': [f"{(count/total_regions)*100:.2f}%" for _, count in remaining_genotypes]
                })
                st.dataframe(df_other_gt, use_container_width=True, hide_index=True)
    
    # Create plots
    col_left, col_right = st.columns(2)
    
    with col_left:
        # Regions by max motif size
        if motif_size_counts:
            # Sort categories logically
            category_order = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10", "Unknown"]
            sorted_counts = {k: motif_size_counts.get(k, 0) for k in category_order if k in motif_size_counts}
            
            df_motif = pd.DataFrame({
                'Max Motif Size': list(sorted_counts.keys()),
                'Number of Regions': list(sorted_counts.values())
            })
            
            # Use distinct theme colors for each bar
            color_palette = get_statistics_color_palette()
            fig_motif = px.bar(
                df_motif, 
                x='Max Motif Size', 
                y='Number of Regions',
                title='Regions by Max Motif Size',
                color='Max Motif Size',
                color_discrete_sequence=color_palette
            )
            fig_motif.update_layout(
                xaxis_title="Max Motif Size (bp)",
                yaxis_title="Number of Regions",
                height=400,
                showlegend=False,
                title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
                xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
                yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
            )
            # Update bar colors to cycle through palette
            for i, trace in enumerate(fig_motif.data):
                trace.marker.color = color_palette[i % len(color_palette)]
            st.plotly_chart(fig_motif, use_container_width=True)
    
    with col_right:
        # Distribution of all motif sizes
        if motif_lengths:
            df_dist = pd.DataFrame({'Motif Size (bp)': motif_lengths})
            color_palette = get_statistics_color_palette()
            fig_dist = px.histogram(
                df_dist,
                x='Motif Size (bp)',
                nbins=20,
                title='Distribution of Motif Sizes',
                labels={'count': 'Frequency'},
                color_discrete_sequence=[color_palette[0]]  # Use primary blue
            )
            fig_dist.update_layout(
                height=400,
                title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
                xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
                yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
            )
            st.plotly_chart(fig_dist, use_container_width=True)
        else:
            st.info("No motif size data available for distribution plot.")
    
    # Regions by chromosome
    if regions_by_chromosome:
        # Sort chromosomes naturally
        sorted_chroms = sorted(
            regions_by_chromosome.keys(),
            key=lambda x: (int(x[3:]) if x[3:].isdigit() else (23 if x == 'chrX' else 24 if x == 'chrY' else 99))
        )
        chrom_data = pd.DataFrame({
            'Chromosome': sorted_chroms,
            'Number of Regions': [regions_by_chromosome[ch] for ch in sorted_chroms]
        })
        
        color_palette = get_statistics_color_palette()
        fig_chrom = px.bar(
            chrom_data,
            x='Chromosome',
            y='Number of Regions',
            title='Regions by Chromosome',
            color='Chromosome',
            color_discrete_sequence=color_palette
        )
        fig_chrom.update_layout(
            xaxis_tickangle=-45,
            height=400,
            showlegend=False,
            title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
            xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
            yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
            xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
            yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
        )
        # Update bar colors to cycle through palette
        for i, trace in enumerate(fig_chrom.data):
            trace.marker.color = color_palette[i % len(color_palette)]
        st.plotly_chart(fig_chrom, use_container_width=True)
    
    # Genotype distribution plot
    if genotype_counts:
        sorted_genotypes = sorted(genotype_counts.items(), 
                                 key=lambda x: (len(x[0].split('/')), x[0]))
        df_genotype = pd.DataFrame({
            'Genotype': [gt for gt, _ in sorted_genotypes],
            'Count': [count for _, count in sorted_genotypes]
        })
        
        color_palette = get_statistics_color_palette()
        fig_genotype = px.bar(
            df_genotype,
            x='Genotype',
            y='Count',
            title='Genotype Distribution',
            color='Genotype',
            color_discrete_sequence=color_palette
        )
        fig_genotype.update_layout(
            xaxis_title="Genotype",
            yaxis_title="Number of Regions",
            height=400,
            showlegend=False,
            title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
            xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
            yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
            xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
            yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
        )
        # Update bar colors to cycle through palette
        for i, trace in enumerate(fig_genotype.data):
            trace.marker.color = color_palette[i % len(color_palette)]
        st.plotly_chart(fig_genotype, use_container_width=True)
    
    st.markdown("---")


@st.cache_data
def compute_sample_statistics(_vcf_file, sample_name):
    """Compute statistics for a single sample VCF file (cached)"""
    stats = {
        'sample_name': sample_name,
        'total_regions': 0,
        'motif_size_counts': {},
        'motif_lengths': [],
        'regions_by_chromosome': {},
        'genotype_counts': {}
    }
    
    try:
        for rec in _vcf_file.fetch():
            stats['total_regions'] += 1
            
            # Get motifs
            motifs = rec.info.get('MOTIFS', [])
            if isinstance(motifs, tuple):
                motifs = list(motifs)
            elif not isinstance(motifs, list):
                motifs = [motifs] if motifs else []
            
            # Calculate max motif size
            max_motif_size = 0
            if motifs:
                motif_lengths_in_region = [len(str(m)) for m in motifs if m]
                if motif_lengths_in_region:
                    max_motif_size = max(motif_lengths_in_region)
                    stats['motif_lengths'].extend(motif_lengths_in_region)
            
            # Categorize
            if max_motif_size == 0:
                category = "Unknown"
            elif max_motif_size <= 10:
                category = str(max_motif_size)
            else:
                category = ">10"
            
            stats['motif_size_counts'][category] = stats['motif_size_counts'].get(category, 0) + 1
            
            # Count by chromosome
            chrom = rec.chrom
            stats['regions_by_chromosome'][chrom] = stats['regions_by_chromosome'].get(chrom, 0) + 1
            
            # Extract genotype information
            if len(rec.samples) > 0 and 'GT' in rec.samples[0]:
                gt = rec.samples[0]['GT']
                if gt is not None:
                    if isinstance(gt, tuple):
                        gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                    elif isinstance(gt, (list,)):
                        gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                    else:
                        gt_str = str(gt)
                    stats['genotype_counts'][gt_str] = stats['genotype_counts'].get(gt_str, 0) + 1
    except Exception as e:
        stats['error'] = str(e)
    
    return stats


@st.cache_data
def compute_sample_statistics_from_path(vcf_file_path, sample_name):
    """Compute statistics for a single sample VCF file using file path (cached)"""
    from collections import Counter
    
    if not vcf_file_path or not Path(vcf_file_path).exists():
        return {'error': 'VCF file not found', 'sample_name': sample_name}
    
    try:
        vcf = pysam.VariantFile(vcf_file_path)
        
        stats = {
            'sample_name': sample_name,
            'total_regions': 0,
            'motif_size_counts': Counter(),
            'motif_lengths': [],
            'regions_by_chromosome': Counter(),
            'genotype_counts': Counter()
        }
        
        for rec in vcf:
            stats['total_regions'] += 1
            
            # Get motifs
            motifs = rec.info.get('MOTIFS', [])
            if isinstance(motifs, tuple):
                motifs = list(motifs)
            elif not isinstance(motifs, list):
                motifs = [motifs] if motifs else []
            
            # Calculate max motif size
            max_motif_size = 0
            if motifs:
                motif_lengths_in_region = [len(str(m)) for m in motifs if m]
                if motif_lengths_in_region:
                    max_motif_size = max(motif_lengths_in_region)
                    stats['motif_lengths'].extend(motif_lengths_in_region)
            
            # Categorize
            if max_motif_size == 0:
                category = "Unknown"
            elif max_motif_size <= 10:
                category = str(max_motif_size)
            else:
                category = ">10"
            
            stats['motif_size_counts'][category] += 1
            
            # Count by chromosome
            chrom = rec.chrom
            stats['regions_by_chromosome'][chrom] += 1
            
            # Extract genotype information
            if len(rec.samples) > 0 and 'GT' in rec.samples[0]:
                gt = rec.samples[0]['GT']
                if gt is not None:
                    if isinstance(gt, tuple):
                        gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                    elif isinstance(gt, (list,)):
                        gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                    else:
                        gt_str = str(gt)
                    stats['genotype_counts'][gt_str] += 1
        
        vcf.close()
        
        return {
            'sample_name': sample_name,
            'total_regions': stats['total_regions'],
            'motif_size_counts': dict(stats['motif_size_counts']),
            'motif_lengths': stats['motif_lengths'],
            'regions_by_chromosome': dict(stats['regions_by_chromosome']),
            'genotype_counts': dict(stats['genotype_counts'])
        }
    except Exception as e:
        return {'error': str(e), 'sample_name': sample_name}


def display_cohort_statistics(cohort_files, cohort_file_paths):
    """Display statistics for cohort VCF files, organized by sample"""
    if not cohort_files or not cohort_file_paths:
        st.info("Please load cohort VCF files to view statistics.")
        return
    
    st.markdown(
        """
        <div style="border-radius: 0.5rem; border: 2px solid #667eea; background: #f7faff; padding: 1.2rem 1.5rem; margin-bottom:0.5rem;">
            <span style="font-size: 1.8rem; font-weight:700; color: #243c5a;">📊 Cohort VCF Statistics (By Sample)</span>
        </div>
        """,
        unsafe_allow_html=True
    )
    
    try:
        # Extract sample names
        sample_names = [f.split(".")[0] if "." in f else f for f in cohort_file_paths]
        
        if not sample_names:
            st.warning("No samples found in cohort files.")
            return
        
        # Initialize session state for selected sample
        if 'selected_sample_idx' not in st.session_state:
            st.session_state.selected_sample_idx = 0
        
        # Sample selection with a nice dropdown
        st.markdown("### 👤 Select Sample to View Statistics")
        
        # Add custom styling for the selectbox
        st.markdown("""
            <style>
            div[data-baseweb="select"] > div {
                font-size: 1.2rem !important;
                padding: 0.75rem !important;
                border-radius: 0.5rem !important;
            }
            div[data-baseweb="select"] > div > div {
                font-size: 1.2rem !important;
                font-weight: 600 !important;
            }
            </style>
        """, unsafe_allow_html=True)
        
        # Use selectbox for sample selection - much cleaner and scrollable
        num_samples = len(sample_names)
        selected_idx = st.selectbox(
            "Choose a sample:",
            options=list(range(num_samples)),
            format_func=lambda x: f"{sample_names[x]} (Sample {x + 1}/{num_samples})",
            index=st.session_state.selected_sample_idx,
            key="cohort_sample_selectbox"
        )
        
        # Update session state
        if selected_idx != st.session_state.selected_sample_idx:
            st.session_state.selected_sample_idx = selected_idx
            st.rerun()
        
        # Get selected sample statistics
        if selected_idx >= len(sample_names):
            selected_idx = 0
            st.session_state.selected_sample_idx = 0
        
        selected_sample_name = sample_names[selected_idx]
        
        # Get base path from session state and construct full path for selected sample
        cohort_base_path = st.session_state.get('path_to_cohort', '')
        if cohort_base_path and not cohort_base_path.endswith('/'):
            cohort_base_path += '/'
        selected_vcf_file_path = cohort_base_path + cohort_file_paths[selected_idx]
        
        # Compute statistics for selected sample only (cached)
        stats = compute_sample_statistics_from_path(selected_vcf_file_path, selected_sample_name)
        
        if stats.get('error'):
            st.error(f"Error computing statistics for {selected_sample_name}: {stats['error']}")
            return
        
        if stats['total_regions'] == 0:
            st.warning(f"No regions found for sample: {selected_sample_name}")
            return
        
        # Display selected sample name prominently
        st.markdown(f"""
            <div style="border-radius: 0.5rem; border: 3px solid #764ba2; background: linear-gradient(135deg, #faf6ff 0%, #f0ebff 100%); padding: 1.5rem; margin: 1rem 0; text-align: center;">
                <span style="font-size: 2rem; font-weight: 800; color: #764ba2;">📊 {selected_sample_name}</span>
            </div>
        """, unsafe_allow_html=True)
        
        # Summary metrics for selected sample
        st.markdown("""
            <style>
            [data-testid="stMetricValue"] {
                font-size: 2rem !important;
            }
            [data-testid="stMetricLabel"] {
                font-size: 1.2rem !important;
            }
            </style>
        """, unsafe_allow_html=True)
        
        total_regions = stats['total_regions']
        motif_size_counts = stats['motif_size_counts']
        motif_lengths = stats['motif_lengths']
        regions_by_chromosome = stats['regions_by_chromosome']
        genotype_counts = stats['genotype_counts']
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Regions", f"{total_regions:,}")
        with col2:
            st.metric("Chromosomes", len(regions_by_chromosome))
        with col3:
            if motif_lengths:
                avg_motif_size = sum(motif_lengths) / len(motif_lengths)
                st.metric("Avg Motif Size", f"{avg_motif_size:.1f} bp")
            else:
                st.metric("Avg Motif Size", "N/A")
        with col4:
            if motif_lengths:
                max_overall = max(motif_lengths)
                st.metric("Max Motif Size", f"{max_overall} bp")
            else:
                st.metric("Max Motif Size", "N/A")
        
        # Genotype statistics
        if genotype_counts:
            with st.container():
                st.markdown(
                    """
                    <div style="border-radius: 0.5rem; border: 2px solid #764ba2; background: #faf6ff; padding: 1.2rem 1.5rem; margin-bottom:0.5rem;">
                        <span style="font-size: 1.75rem; font-weight: 700; color: #4B247B;">🧬 Genotype Distribution</span>
                    </div>
                    """,
                    unsafe_allow_html=True
                )
            sorted_genotypes = sorted(genotype_counts.items(), 
                                     key=lambda x: (len(x[0].split('/')), x[0]))
            
            num_genotype_cols = min(len(sorted_genotypes), 6)
            if num_genotype_cols > 0:
                genotype_cols = st.columns(num_genotype_cols)
                for idx, (gt, count) in enumerate(sorted_genotypes[:num_genotype_cols]):
                    with genotype_cols[idx]:
                        percentage = (count / total_regions) * 100
                        st.metric(f"Genotype {gt}", f"{count:,}", f"{percentage:.1f}%")
                
                if len(sorted_genotypes) > 6:
                    remaining_genotypes = sorted_genotypes[6:]
                    st.markdown("""
                        <div style="font-size: 1.25rem; font-weight: 600; margin-top: 1rem; margin-bottom: 0.5rem;">
                            <strong>Other genotypes:</strong>
                        </div>
                    """, unsafe_allow_html=True)
                    df_other_gt = pd.DataFrame({
                        'Genotype': [gt for gt, _ in remaining_genotypes],
                        'Count': [count for _, count in remaining_genotypes],
                        'Percentage': [f"{(count/total_regions)*100:.2f}%" for _, count in remaining_genotypes]
                    })
                    st.dataframe(df_other_gt, use_container_width=True, hide_index=True)
        
        # Create plots for selected sample
        col_left, col_right = st.columns(2)
        
        with col_left:
            if motif_size_counts:
                category_order = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10", "Unknown"]
                sorted_counts = {k: motif_size_counts.get(k, 0) for k in category_order if k in motif_size_counts}
                
                df_motif = pd.DataFrame({
                    'Max Motif Size': list(sorted_counts.keys()),
                    'Number of Regions': list(sorted_counts.values())
                })
                
                color_palette = get_statistics_color_palette()
                fig_motif = px.bar(
                    df_motif,
                    x='Max Motif Size',
                    y='Number of Regions',
                    title=f'Regions by Max Motif Size - {selected_sample_name}',
                    color='Max Motif Size',
                    color_discrete_sequence=color_palette
                )
                fig_motif.update_layout(
                    xaxis_title="Max Motif Size (bp)",
                    yaxis_title="Number of Regions",
                    height=400,
                    showlegend=False,
                    title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
                    xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                    yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                    xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
                    yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
                )
                for i, trace in enumerate(fig_motif.data):
                    trace.marker.color = color_palette[i % len(color_palette)]
                st.plotly_chart(fig_motif, use_container_width=True)
        
        with col_right:
            if motif_lengths:
                df_dist = pd.DataFrame({'Motif Size (bp)': motif_lengths})
                color_palette = get_statistics_color_palette()
                fig_dist = px.histogram(
                    df_dist,
                    x='Motif Size (bp)',
                    nbins=20,
                    title=f'Distribution of Motif Sizes - {selected_sample_name}',
                    labels={'count': 'Frequency'},
                    color_discrete_sequence=[color_palette[0]]
                )
                fig_dist.update_layout(
                    height=400,
                    title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
                    xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                    yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                    xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
                    yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
                )
                st.plotly_chart(fig_dist, use_container_width=True)
        
        # Regions by chromosome
        if regions_by_chromosome:
            sorted_chroms = sorted(
                regions_by_chromosome.keys(),
                key=lambda x: (int(x[3:]) if x[3:].isdigit() else (23 if x == 'chrX' else 24 if x == 'chrY' else 99))
            )
            chrom_data = pd.DataFrame({
                'Chromosome': sorted_chroms,
                'Number of Regions': [regions_by_chromosome[ch] for ch in sorted_chroms]
            })
            
            color_palette = get_statistics_color_palette()
            fig_chrom = px.bar(
                chrom_data,
                x='Chromosome',
                y='Number of Regions',
                title=f'Regions by Chromosome - {selected_sample_name}',
                color='Chromosome',
                color_discrete_sequence=color_palette
            )
            fig_chrom.update_layout(
                xaxis_tickangle=-45,
                height=400,
                showlegend=False,
                title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
                xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
                yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
            )
            for i, trace in enumerate(fig_chrom.data):
                trace.marker.color = color_palette[i % len(color_palette)]
            st.plotly_chart(fig_chrom, use_container_width=True)
        
        # Genotype distribution plot
        if genotype_counts:
            sorted_genotypes = sorted(genotype_counts.items(), 
                                     key=lambda x: (len(x[0].split('/')), x[0]))
            df_genotype = pd.DataFrame({
                'Genotype': [gt for gt, _ in sorted_genotypes],
                'Count': [count for _, count in sorted_genotypes]
            })
            
            color_palette = get_statistics_color_palette()
            fig_genotype = px.bar(
                df_genotype,
                x='Genotype',
                y='Count',
                title=f'Genotype Distribution - {selected_sample_name}',
                color='Genotype',
                color_discrete_sequence=color_palette
            )
            fig_genotype.update_layout(
                xaxis_title="Genotype",
                yaxis_title="Number of Regions",
                height=400,
                showlegend=False,
                title_font=dict(size=26, color='black', family='Arial Black, sans-serif'),
                xaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                yaxis_title_font=dict(size=22, color='black', family='Arial Black, sans-serif'),
                xaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif'),
                yaxis_tickfont=dict(size=20, color='black', family='Arial Black, sans-serif')
            )
            for i, trace in enumerate(fig_genotype.data):
                trace.marker.color = color_palette[i % len(color_palette)]
            st.plotly_chart(fig_genotype, use_container_width=True)
        
        st.markdown("---")
        
    except Exception as e:
        st.error(f"Error computing cohort statistics: {e}")


def render_sidebar_branding():

    try:
        
        logo_path = Path(__file__).resolve().parent.parent / "ProleTRact_logo.svg"
        if logo_path.exists():
            # Read the image as base64 for browser embedding
            import base64
            with open(logo_path, "rb") as image_file:
                encoded = base64.b64encode(image_file.read()).decode("utf-8")
            st.sidebar.markdown(
                f"""
                <div style="display: flex; justify-content: center; align-items: center; margin-bottom: 2.2em;">
                    <img src="data:image/svg+xml;base64,{encoded}" alt="ProleTRact Logo" style="max-width:250px; width:250px; display:block;" />
                </div>
                """,
                unsafe_allow_html=True
            )
    except FileNotFoundError:
        st.write("Logo not found")



    # Page links directly in the sidebar
    with st.sidebar:
        # Add extra padding and make the font size *twice* as big for sidebar page links
        st.markdown("""
            <style>
                .stSidebarNav .stPageLink {
                    padding: 1.8rem 2.5rem 1.8rem 2.5rem !important;
                    font-size: 4.2rem !important;
                    border-radius: 15px;
                    min-height: 92px;
                    margin-bottom: 9px;
                }
                .stSidebarNav .stPageLink:hover {
                    background: #ece9ff !important;
                }
            </style>
        """, unsafe_allow_html=True)

        st.page_link("app.py", label="About ProleTRact", icon="🏠")
        st.page_link("pages/1_Individual_sample.py", label="Individual Mode 👤")
        st.page_link("pages/2_Cohort_Reads.py", label="Cohort Mode (Reads-based)    ☰")
        st.page_link("pages/3_Cohort_Assembly.py", label="Cohort Mode (Assembly-based)   ━━━━") #━━━━


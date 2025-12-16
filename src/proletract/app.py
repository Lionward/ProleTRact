import streamlit as st
import pandas as pd
import io
import numpy as np
import altair as alt
import sys
from pathlib import Path

try:
    from proletract.app_shared import configure_page, ensure_state_initialized
    from proletract.modules.viz.vis_helper import (
        create_genotype_comparison_matrix,
        display_dynamic_sequence_with_highlighted_motifs,
        motif_legend_html,
    )
    from proletract.modules.viz import utils
    from proletract.modules.viz import plots
except ModuleNotFoundError:
    src_path = Path(__file__).resolve().parents[1]
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))

    from proletract.app_shared import configure_page, ensure_state_initialized
    from proletract.modules.viz.vis_helper import (
        create_genotype_comparison_matrix,
        display_dynamic_sequence_with_highlighted_motifs,
        motif_legend_html,
    )
    from proletract.modules.viz import utils
    from proletract.modules.viz import plots


if __package__ is None or __package__ == "":
    src_path = Path(__file__).resolve().parents[1]
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))


def main():
    configure_page()
    ensure_state_initialized()
    st.subheader("ProleTRact: Tandem Repeat Analysis Portal")
    st.info("🔍 Explore, visualize, and compare tandem repeat regions from TandemTwister outputs. Use the pages on the sidebar to start with an individual sample or cohort.")

    quickstart_tab, examples_tab, faq_tab, upcoming_tab = st.tabs(["Quickstart", "Examples", "FAQ", "Upcoming Features"])

    with quickstart_tab:
        st.markdown("### Get started in 3 steps")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("**1) Individual**")
            st.caption("Load one VCF and explore TR regions interactively.")
            st.page_link("pages/1_Individual_sample.py", label="Open Individual page", icon="👤")
        with col2:
            st.markdown("**2) Cohort (Reads)**")
            st.caption("Aggregate read-based VCFs and compare across samples.")
            st.page_link("pages/2_Cohort_Reads.py", label="Open Cohort Mode (Reads-based) 👤👤👤")
        with col3:
            st.markdown("**3) Cohort (Assembly)**")
            st.caption("Analyze haplotype-resolved assembly VCFs.")
            st.page_link("pages/3_Cohort_Assembly.py", label="Open Cohort Mode (Assembly-based) 👤👤👤")

        st.markdown("---")
       

    with examples_tab:


 
        st.markdown("---")
        st.markdown("### Plot gallery (same components used in the tool)")

        # A) Sequence with highlighted motifs
        st.markdown("**A) Sequence with highlighted motifs**")
        demo_motif_names = ["CAG", "GAA", "TTC"]
        demo_colors = {i: c for i, c in enumerate(utils.get_color_palette(len(demo_motif_names)))}
        # Build a synthetic sequence with segments: CAGx10, GAAx6, TTCx8 (interruptions inserted)
        seq = "CAG" * 10 + "A" + "GAA" * 6 + "TT" + "TTC" * 8
        # spans are per motif occurrence: each individual motif gets its own span
  
        spans = "(1-3),(4-6),(7-9),(10-12),(13-15),(16-18),(19-21),(22-24),(25-27),(28-30),(32-34),(35-37),(38-40),(41-43),(44-46),(47-49),(52-54),(55-57),(58-60),(61-63),(64-66),(67-69),(70-72),(73-75)"
        motif_ids = [0] * 10 + [1] * 6 + [2] * 8  # 10 CAG, 6 GAA, 8 TTC
           # Show the motif legend using the same function
        motif_legend_html(motif_ids, demo_colors, demo_motif_names)
        st.caption(
            "Legend explains which colors correspond to which motifs in the region. Also displays summary statistics."
        )

        st.markdown("---")
        display_dynamic_sequence_with_highlighted_motifs(
            "Demo", seq, motif_ids, spans, demo_colors, demo_motif_names, supporting_reads=None
        )
        st.caption(
            "Each colored block corresponds to a detected motif run; interruptions (red highlighted nucleotides) are shown between runs."
            " This mirrors the individual-sample sequence panel."
        )

        st.markdown("**Motif legend (used across visualizations)**")
     
        # B) Genotype comparison matrix
        st.markdown("**B) Genotype comparison matrix**")
        demo_genotypes = {
            "case_01": "0/1",
            "case_02": "1/1",
            "ctrl_01": "0/0",
            "ctrl_02": "0/1",
            "ctrl_03": "0/0",
        }
        create_genotype_comparison_matrix(demo_genotypes)
        st.caption("Summarizes genotypes across samples with color- and icon-coding for quick scanning.")

        st.markdown("---")
        # C) Stack plot + heatmap
        st.markdown("**C) Stack plot and heatmap of motif segments**")
        # Build minimal structures expected by stack_plot
        record = {
            'motifs': demo_motif_names,
            'chr': 'chr1',
            'pos': 1000,
            'stop': 3000,
            'id': 'chr1:1000-3000'
        }
        sequences = [
            {'name': 'S1_alle1', 'sequence': "CAG"*8 + "GAA"*3 + "TTC"*2},
            {'name': 'S1_alle2', 'sequence': "CAG"*12 + "TTC"*2},
            {'name': 'S2_alle1', 'sequence': "GAA"*5 + "CAG"*6},
            {'name': 'S3_alle1', 'sequence': "TTC"*40 + "GAA"*30 + "CAG"*20},
        ]
        span_list = [
            "(1-24),(25-33),(34-39)",
            "(1-36),(37-42)",
            "(1-15),(16-33)",
            "(1-120),(121-210),(211-270)",
        ]
        motif_ids_list = [
            [0, 1, 2],
            [0, 2],
            [1, 0],
            [2, 1, 0],
        ]
        # Render the stack plot + internal heatmap and summary stats
        _motif_colors, df_stack = plots.stack_plot(record, demo_motif_names, sequences, span_list, motif_ids_list, sort_by="Value", max_height=600, max_width=800)
        st.caption(
            "Stack plot: each row is a sample/allele; colored blocks represent motif runs along the sequence."
            " The heatmap on top aggregates motif occurrences by sample and motif."
        )

        st.markdown("---")
        # D) Motif count per sample (cohort bar chart)
        st.markdown("**D) Motif count per sample (with threshold)**")
        # Build a minimal dataframe 
        motif_names = demo_motif_names
        samples = ["S1_alle1", "S1_alle2", "S2_alle1", "S2_alle2", "S3_alle1"]
        rows = []
        rng = np.random.default_rng(7)
        for s in samples:
            n = int(rng.integers(5, 13))
            for _ in range(n):
                rows.append({"Sample": s, "Motif": rng.choice(motif_names)})
        demo_df = pd.DataFrame(rows)
        region = st.session_state.pathogenic_TRs.iloc[0]['region'] if "pathogenic_TRs" in st.session_state and not st.session_state.pathogenic_TRs.empty else "chr1:1000-2000"
        plots.bar_plot_motif_count(demo_df, region, sort_by="Value")
        st.caption(
            "Bar chart: total motif segments per sample/allele. If the displayed region is in the pathogenic catalog, a red threshold line is shown."
        )
        st.markdown("---")
        st.markdown("### Example of pathogenic TRs input file")
        if "pathogenic_TRs" in st.session_state:
            df = st.session_state.pathogenic_TRs
            preview = df[["chrom", "start", "end", "motif", "disease", "gene"]].head(10)
            # Create the custom HTML table with larger font sizes
            table_html = """
            <style>
                .pathogenic-table {
                    width: 100%;
                    border-collapse: collapse;
                    margin: 20px 0;
                    font-size: 24px !important;
                }
                .pathogenic-table th {
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                    padding: 16px 12px;
                    text-align: left;
                    font-size: 26px !important;
                    font-weight: bold !important;
                    border: 1px solid #e2e8f0;
                }
                .pathogenic-table td {
                    padding: 14px 12px;
                    border: 1px solid #e2e8f0;
                    font-size: 24px !important;
                    background-color: #f9fafb;
                }
                .pathogenic-table tr:nth-child(even) td {
                    background-color: #ffffff;
                }
                .pathogenic-table tr:hover td {
                    background-color: #f1f5f9;
                }
            </style>
            <table class="pathogenic-table">
                <thead>
                    <tr>
            """
            # Add the header row
            for col in preview.columns:
                table_html += f"<th>{col}</th>"
            table_html += """
                    </tr>
                </thead>
                <tbody>
            """
            # Add the data rows
            for _, row in preview.iterrows():
                table_html += "<tr>"
                for col in preview.columns:
                    table_html += f"<td>{row[col]}</td>"
                table_html += "</tr>"
            table_html += """
                </tbody>
            </table>
            """
            st.markdown(table_html, unsafe_allow_html=True)
        else:
            st.warning("Pathogenic TR catalog not loaded.")

      
    with upcoming_tab:
        st.markdown("### 🚀 Upcoming Features")
        st.markdown("We're continuously improving ProleTRact. Here are some features planned for future releases:")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### 🔬 Read-level visualization")
            st.info(
                "**IGV-like read browser**: Interactive read alignment viewer showing individual reads mapped to TR regions. "
                "Visualize read depth, support for different alleles, and identify read-level evidence for expansions/contractions. "
                "Zoom, pan, and filter reads by quality or allele support.",
                icon="📊"
            )
            
            st.markdown("#### 🔍 Advanced filtering")
            st.info(
                "**Multi-criteria filtering**: Filter samples/regions by repeat count ranges, quality scores, coverage depth, "
                "pathogenicity thresholds",
                icon="🎯"
            )
            
            
            st.markdown("#### 🧬 Phasing & inheritance")
            st.info(
                "**Family pedigree visualization**: Track TR alleles through pedigrees. Visualize segregation patterns, "
                "identify de novo expansions, and analyze inheritance modes. Support for trio and extended family analysis.",
                icon="👨‍👩‍👧‍👦"
            )
            
        
        with col2:
            st.markdown("#### 🌐 Population databases")
            st.info(
                "**Integration with gnomAD & HGSVC**: Compare sample alleles against population reference databases. "
                "Display allele frequencies from healthy cohorts. Identify rare or novel expansions relative to population data.",
                icon="🌍"
            )
            
            st.markdown("#### 🧪 Quality metrics")
            st.info(
                "**Comprehensive QC dashboard**: Visualize per-sample and per-region quality metrics (coverage, depth, allelic balance). "
                "Identify problematic regions or samples. Automated QC flagging and quality score filtering.",
                icon="✅"
            )
            
            st.markdown("#### 📈 Statistical analysis")
            st.info(
                "**Allele frequency analysis**: Calculate population allele frequencies, compare case vs control distributions, "
                "perform statistical tests (t-test, Mann-Whitney, etc.). Generate summary statistics and effect size estimates.",
                icon="📉"
            )
        st.markdown("---")
        st.markdown("#### 💡 Have suggestions?")
        st.markdown(
            "We'd love to hear your ideas! Contact us or open an issue on GitHub to suggest features or report bugs. "
            "Your feedback helps shape the future of ProleTRact."
        )
        

    with faq_tab:
        with st.expander("What input formats are supported?", expanded=False):
            st.info("VCF files produced by TandemTwister (reads-based or assembly-based).")
        with st.expander("How do I select a TR region?", expanded=False):
            st.info("On the Individual page/cohort pages, search by coordinate and then visualize.")
        with st.expander("How do I group samples in a cohort?", expanded=False):
            st.info("Run TandemTwister on a cohort of samples and then upload the VCF files to the cohort pages by specifying the path to the VCF files.")
        with st.expander("How do I use the pathogenic TR catalog?", expanded=False):
            st.info("Run TandemTwister using the pathogenic TR catalog as the input file and then upload the VCF files to the cohort pages by specifying the path to the VCF files.")
if __name__ == "__main__":
    main()
    

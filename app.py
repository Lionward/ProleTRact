import streamlit as st
import pandas as pd
import io
import numpy as np
import altair as alt
from proletract.app_shared import configure_page, ensure_state_initialized
from proletract.modules.viz.visualization import Visualization
from proletract.modules.viz.vis_helper import create_genotype_comparison_matrix, display_dynamic_sequence_with_highlighted_motifs


def main():
    configure_page()
    ensure_state_initialized()

    st.title("🧬 ProleTRact")
    st.subheader("Tandem Repeat Analysis Portal")

    st.info(
        "Explore, visualize, and compare tandem repeat regions from TandemTwister outputs."
        " Use the pages on the sidebar to start with an individual sample or cohort.",
        icon="💡",
    )

    quickstart_tab, examples_tab, faq_tab = st.tabs(["Quickstart", "Examples", "FAQ"])

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
        viz = Visualization()

        # A) Sequence with highlighted motifs (same renderer as individual view)
        st.markdown("**A) Sequence with highlighted motifs**")
        demo_motif_names = ["CAG", "GAA", "TTC"]
        demo_colors = {i: c for i, c in enumerate(viz.get_color_palette(len(demo_motif_names)))}
        # Build a synthetic sequence with segments: CAGx10, GAAx6, TTCx8 (interruptions inserted)
        seq = "CAG" * 10 + "A" + "GAA" * 6 + "TT" + "TTC" * 8
        # spans are 1-based inclusive pairs, matching motif_ids order below
        spans = "(1-30),(32-49),(52-75)"
        motif_ids = [0, 1, 2]  # indices into demo_motif_names
        display_dynamic_sequence_with_highlighted_motifs(
            "Demo", seq, motif_ids, spans, demo_colors, demo_motif_names, supporting_reads=None
        )
        st.caption(
            "Each colored block corresponds to a detected motif run; interruptions (red highlighted nucleotides) are shown between runs."
            " This mirrors the individual-sample sequence panel."
        )

        st.markdown("---")
        # B) Genotype comparison matrix (as used in cohort visualizations)
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
        # C) Stack plot + heatmap (same routine used in cohort stack visualization)
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
            # Make this sequence much bigger by repeating motifs more
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
        # Render stack plot + internal heatmap and summary stats
        _motif_colors, df_stack = viz.stack_plot(record, demo_motif_names, sequences, span_list, motif_ids_list, sort_by="Value")
        st.caption(
            "Stack plot: each row is a sample/allele; colored blocks represent motif runs along the sequence."
            " The heatmap on top aggregates motif occurrences by sample and motif."
        )

        st.markdown("---")
        # D) Motif count per sample (cohort bar chart)
        st.markdown("**D) Motif count per sample (with threshold)**")
        # Build a minimal dataframe expected by bar_plot_motif_count
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
        viz.bar_plot_motif_count(demo_df, region, sort_by="Value")
        st.caption(
            "Bar chart: total motif segments per sample/allele. If the displayed region is in the pathogenic catalog, a red threshold line is shown."
        )
        st.markdown("---")
        st.markdown("### Example of pathogenic TRs input file")
        if "pathogenic_TRs" in st.session_state:
            df = st.session_state.pathogenic_TRs
            preview = df[["chrom", "start", "end", "motif", "disease", "gene"]].head(10)
            st.dataframe(preview, use_container_width=True, hide_index=True)
        else:
            st.warning("Pathogenic TR catalog not loaded.")

      
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
    

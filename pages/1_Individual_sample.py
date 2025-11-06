import streamlit as st
from proletract.app_shared import configure_page, ensure_state_initialized, display_vcf_statistics


def run_page():
    configure_page()
    ensure_state_initialized()

    vcf_handler = st.session_state.vcf_handler
    visualization = st.session_state.visualization

    st.header("👤 Individual sample")
    vcf_handler.handle_individual_sample()
    st.session_state.analysis_mode = "individual sample"

    # Check if VCF is loaded
    vcf_loaded = 'vcf_file_path' in st.session_state and st.session_state.vcf_file_path
    records_loaded = 'records_map' in st.session_state

    # Create tabs for statistics and visualization
    if vcf_loaded or records_loaded:
        stats_tab, viz_tab = st.tabs(["📊 Statistics", "🔬 Region Analysis"])
        
        with stats_tab:
            if vcf_loaded:
                display_vcf_statistics(st.session_state.vcf_file_path, mode='individual')
            else:
                st.info("Please load a VCF file to view statistics.")
        
        with viz_tab:
            if records_loaded:
                visualization.visulize_region()
            else:
                st.info("Please load a VCF file to start region analysis.")
    else:
        st.info(
            "👈 Please load a VCF file from the sidebar to get started.\n\n"
            "Additionally, ensure that your cohort folder contains the same repeat catalog that was used for this sample."
        )
        st.stop()


if __name__ == "__main__":
    run_page()


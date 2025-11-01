import streamlit as st
from proletract.app_shared import configure_page, ensure_state_initialized


def run_page():
    configure_page()
    ensure_state_initialized()

    vcf_handler = st.session_state.vcf_handler
    visualization = st.session_state.visualization

    st.header("👤 Individual sample")
    vcf_handler.handle_individual_sample()
    st.session_state.analysis_mode = "individual sample"

    if 'records_map' in st.session_state:
        visualization.visulize_region()
    else:
        st.stop()


if __name__ == "__main__":
    run_page()


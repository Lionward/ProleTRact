import streamlit as st
from proletract.app_shared import configure_page, ensure_state_initialized


def run_page():
    configure_page()
    ensure_state_initialized()

    cohort_handler = st.session_state.cohort_handler
    visualization = st.session_state.visualization

    st.header("👥 Cohort – Reads-based VCF")
    st.session_state['cohort_mode'] = "reads"
    cohort_handler.handle_cohort()
    st.session_state.analysis_mode = "Cohort"
    visualization.visulize_cohort()


if __name__ == "__main__":
    run_page()


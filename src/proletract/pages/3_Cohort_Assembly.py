import sys
from pathlib import Path

if __package__ is None or __package__ == "":
    src_path = Path(__file__).resolve().parents[2]
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))

import streamlit as st

try:
    from proletract.app_shared import (
        configure_page,
        ensure_state_initialized,
        display_cohort_statistics,
    )
except ModuleNotFoundError:  # pragma: no cover - local dev fallback
    import sys
    from pathlib import Path

    src_path = Path(__file__).resolve().parents[2]
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))

    from proletract.app_shared import (
        configure_page,
        ensure_state_initialized,
        display_cohort_statistics,
    )


def run_page():
    configure_page()
    ensure_state_initialized()

    cohort_handler = st.session_state.cohort_handler
    visualization = st.session_state.visualization

    st.header("👥 Cohort – Assembly VCF")
    st.session_state['cohort_mode'] = "assembly"
    cohort_handler.handle_cohort()
    st.session_state.analysis_mode = "Cohort"
    
    # Check if cohort files are loaded
    cohort_loaded = ('cohort_files' in st.session_state and 'cohort_file_paths' in st.session_state and
                     st.session_state.cohort_files and st.session_state.cohort_file_paths)
    
    # Create tabs for statistics and visualization
    if cohort_loaded:
        stats_tab, viz_tab = st.tabs(["📊 Statistics", "🔬 Cohort Analysis"])
        
        with stats_tab:
            display_cohort_statistics(st.session_state.cohort_files, st.session_state.cohort_file_paths)
        
        with viz_tab:
            visualization.visulize_cohort()
    else:
        st.info("👈 Please load cohort VCF files from the sidebar to get started.")


if __name__ == "__main__":
    run_page()


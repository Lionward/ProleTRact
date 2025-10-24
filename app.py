import streamlit as st
import streamlit.components.v1 as components
from data_handling import VCFHandler, CohortHandler
from visualization import Visualization
import pandas as pd 
import base64
import os



def get_pathogenic_TRs():
    current_path = os.getcwd()
    pathogenic_trs = pd.read_csv(current_path+"/pathogenic_TRs.bed", sep="\t", header=None)
    pathogenic_trs.columns = ["chrom",'start','end','motif','pathogenic_min','inheritance','disease','gene']
    st.session_state.pathogenic_TRs = pathogenic_trs
    pathogenic_trs['region'] = pathogenic_trs['chrom'] + ":" + pathogenic_trs['start'].astype(str) + "-" + pathogenic_trs['end'].astype(str)


def main():

    if 'vcf_handler' not in st.session_state:
        st.session_state.vcf_handler = VCFHandler()
    if 'cohort_handler' not in st.session_state:
        st.session_state.cohort_handler = CohortHandler()
    st.session_state.visualization = Visualization()
    get_pathogenic_TRs()
    
    vcf_handler = st.session_state.vcf_handler
    cohort_handler = st.session_state.cohort_handler
    visualization = st.session_state.visualization

 
        
    analysis_mode = st.sidebar.radio(
        "Select the type of analysis", 
        ("individual sample", "Cohort"), 
        key="analysis_mode_radio",
        label_visibility='visible',
        help="Choose the analysis workflow.",
    )
    st.markdown("""
        <style>
            /* Make radiogroup container and background fully transparent */
            [data-testid="stSidebar"] div[role="radiogroup"] {
                background: transparent !important;
                box-shadow: none !important;
                border: none !important;
                padding: 0 !important;
                margin: 0 !important;
            }
            /* Make internal radio controls and labels transparent/native and styled */
            [data-testid="stSidebar"] div[data-baseweb="radio"] {
                background: transparent !important;
                box-shadow: none !important;
            }
            [data-testid="stSidebar"] div[data-baseweb="radio"] label {
                font-weight: 600 !important;
                color: #374151 !important;
                background: transparent !important;
            }
            [data-testid="stSidebar"] div[data-baseweb="radio"] > div {
                margin-bottom: 10px !important;
                background: transparent !important;
            }
            /* Remove any extra container backgrounds */
            [data-testid="stSidebar"] div[data-testid='stMarkdownContainer'] {
                background: transparent !important;
                border-radius: 10px;
                padding: 2px 5px;
                margin-bottom: 2px;
            }
        </style>
    """, unsafe_allow_html=True)

    if analysis_mode == "individual sample":
        vcf_handler.handle_individual_sample()
        st.session_state.analysis_mode = "individual sample"
        
        if 'records_map' in st.session_state:
            
            visualization.visulize_region() 
        else:
            st.stop()


    elif analysis_mode == "Cohort":
        cohort_handler.handle_cohort()
        st.session_state.analysis_mode = "Cohort"
        # if the path doesn't end with a slash add it

        visualization.visulize_cohort()


if __name__ == "__main__":

    
    st.set_page_config(layout="wide")
    # Separate button size styling for sidebar vs. main area
    st.markdown("""
        <style>
            /* Main area buttons (make smaller) */
            div[data-testid="stAppViewContainer"] .stButton > button {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                border: none;
                padding: 4px 10px; /* Smaller vertical & horizontal padding */
                border-radius: 10px;
                font-weight: 700 !important;
                font-size: 13px !important;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
                margin: 4px;
                min-width: 110px;
                min-height: 38px;
                letter-spacing: 0.65px;
                line-height: 1.15 !important;
            }
            div[data-testid="stAppViewContainer"] .stButton > button * {
                font-size: 16px !important;
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
            /* Sidebar buttons (still smaller, compact) */
            section[data-testid="stSidebar"] .stButton > button {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                border: none;
                padding: 4px 12px;
                border-radius: 10px;
                font-weight: 600 !important;
                font-size: 13px !important;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                box-shadow: 0 2px 6px rgba(102, 126, 234, 0.12);
                margin: 3px 0;
                min-width: 80px;
                min-height: 28px;
                letter-spacing: 0.3px;
                line-height: 1.1 !important;
            }
            section[data-testid="stSidebar"] .stButton > button * {
                font-size: 13px !important;
                font-weight: 600 !important;
                line-height: 1.1 !important;
            }
            section[data-testid="stSidebar"] .stButton > button:hover {
                transform: translateY(-1.5px) scale(1.01);
                box-shadow: 0 3px 8px rgba(102, 126, 234, 0.15);
                background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
                font-size: 14px !important;
            }
            section[data-testid="stSidebar"] .stButton > button:hover * {
                font-size: 14px !important;
            }
            section[data-testid="stSidebar"] .stButton > button:active,
            section[data-testid="stSidebar"] .stButton > button:focus {
                transform: scale(0.97);
                box-shadow: 0 1px 3px rgba(102, 126, 234, 0.10);
                font-size: 13px !important;
            }
            section[data-testid="stSidebar"] .stButton > button:active *,
            section[data-testid="stSidebar"] .stButton > button:focus * {
                font-size: 13px !important;
            }
            section[data-testid="stSidebar"] .stButton button > div, 
            section[data-testid="stSidebar"] .stButton button > span, 
            section[data-testid="stSidebar"] .stButton button > div > span {
                font-size: 13px !important;
                font-weight: 600 !important;
                line-height: 1.1 !important;
            }
        </style>
    """, unsafe_allow_html=True)
    
    # Add JavaScript to make navigation buttons fixed
    components.html("""
        <script>
        (function() {
            function makeNavButtonsFixed() {
                const buttons = document.querySelectorAll('button');
                let prevButton = null;
                let nextButton = null;
                
                buttons.forEach(btn => {
                    const text = btn.textContent || btn.innerText || '';
                    if (text.includes('Previous region') && btn.style.position !== 'fixed') {
                        prevButton = btn;
                    }
                    if (text.includes('Next region') && btn.style.position !== 'fixed') {
                        nextButton = btn;
                    }
                });
                
                if (prevButton && nextButton) {
                    // Style the buttons to be fixed
                    prevButton.style.position = 'fixed';
                    prevButton.style.bottom = '20px';
                    prevButton.style.left = 'calc(50% - 150px)';
                    prevButton.style.zIndex = '999';
                    
                    nextButton.style.position = 'fixed';
                    nextButton.style.bottom = '20px';
                    nextButton.style.left = 'calc(50% + 50px)';
                    nextButton.style.zIndex = '999';
                }
            }
            
            // Run immediately
            makeNavButtonsFixed();
            
            // Run after a delay
            setTimeout(makeNavButtonsFixed, 100);
            setTimeout(makeNavButtonsFixed, 500);
            setTimeout(makeNavButtonsFixed, 1000);
            
            // Use MutationObserver to watch for DOM changes
            const observer = new MutationObserver(function(mutations) {
                makeNavButtonsFixed();
            });
            
            observer.observe(document.body, {
                childList: true,
                subtree: true
            });
        })();
        </script>
    """, height=0)
    
    st.sidebar.markdown(
        """
        <div style="
            text-align: center;
            font-size: 2.3rem;
            font-weight: 800;
            letter-spacing: 1px;
            margin-bottom: -0.3em;
        "><span style="color:white; font-weight:bold;">🧬 ProleTRact</span></div>
        <div style="color:white; text-align:center; font-size:14px; margin-top:-8px; margin-bottom: 12px;">
            Tandem Repeat Analysis Portal
        </div>
        """, unsafe_allow_html=True
    )
    with st.sidebar.expander("✨ About ProleTRact", expanded=True):
        st.info(
            "Welcome! Explore, visualize, and compare tandem repeat data in samples and cohorts. Switch modes below to begin. 🌀",
            icon="💡"
        )


    main()
    

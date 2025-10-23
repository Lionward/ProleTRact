import streamlit as st
from data_handling import VCFHandler, CohortHandler
from visualization import Visualization
import pandas as pd 
import base64
import os

def get_image_base64(path):

    with open(path, "rb") as f:
        encoded_string = base64.b64encode(f.read()).decode()
    return encoded_string

def get_pathogenic_TRs():
    current_path = os.getcwd()
    pathogenic_trs = pd.read_csv(current_path+"/pathogenic_TRs.bed", sep="\t", header=None)
    pathogenic_trs.columns = ["chrom",'start','end','motif','pathogenic_min','inheritance','disease','gene']
    st.session_state.pathogenic_TRs = pathogenic_trs
    # make column region
    pathogenic_trs['region'] = pathogenic_trs['chrom'] + ":" + pathogenic_trs['start'].astype(str) + "-" + pathogenic_trs['end'].astype(str)


def main():

    if 'vcf_handler' not in st.session_state:
        st.session_state.vcf_handler = VCFHandler()
    if 'cohort_handler' not in st.session_state:
        st.session_state.cohort_handler = CohortHandler()
    #if 'visualization' not in st.session_state:
    st.session_state.visualization = Visualization()
    get_pathogenic_TRs()
    
    vcf_handler = st.session_state.vcf_handler
    cohort_handler = st.session_state.cohort_handler
    visualization = st.session_state.visualization

 
        
    analysis_mode = st.sidebar.radio(
        "Select the type of analysis", 
        ("individual sample", "Cohort", "comparison"), 
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
        # posistion the button in the center
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
    elif analysis_mode == "comparison":
        if 'all_files_parsed' not in st.session_state:
            st.session_state.all_files_parsed = False
        st.session_state.analysis_mode = "comparison"
        vcf_handler.handle_comparison_samples()
        if st.session_state.all_files_parsed:
            visualization.compare_different_technologies()
        else:
            st.warning("Failed to parse the VCF files")
            #st.stop()  

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

    # Path to your local logo
    logo_path = "tandem_twister_vis_logo.png"
    logo_base64 = get_image_base64(logo_path)
    # Custom HTML/CSS to align the logo to the left with a transparent background
    st.sidebar.markdown(
        f"""
        <style>
        .sidebar {{
            background: linear-gradient(180deg, #667eea 0%, #764ba2 100%);
            padding: 0;
            margin: 0;
        }}
        .sidebar .sidebar-content {{
            background: transparent !important;
        }}
        .proletract-title-container {{
            text-align: center;
            margin-top: 32px;
            margin-bottom: -8px;
        }}
        .proletract-title {{
            font-size: 2.3rem;
            font-weight: 800;
            letter-spacing: 1px;
            color: white;
            font-family: inherit;
            display: inline-block;
            line-height: 1.12;
        }}
        .proletract-subtitle {{
            color: white;
            text-align: center;
            font-size: 14px;
            margin-top: -8px;
            margin-bottom: 10px;
            font-family: inherit;
            font-weight: 500;
        }}
        .logo-container {{
            display: flex;
            justify-content: center;
            align-items: center;
            padding: 25px 20px;
            background: rgba(255, 255, 255, 0.0); /* Make logo container fully transparent */
            backdrop-filter: none;
            border-bottom: 1px solid rgba(255, 255, 255, 0.2);
            margin-bottom: 25px;
        }}
        .logo {{
            display: flex;
            justify-content: center;
            align-items: center;
            background: transparent; /* Remove white frame and make background transparent */
            padding: 0;              /* Remove extra padding so image stands alone */
            border-radius: 0;        /* Remove border radius so image outline is unchanged */
            box-shadow: none;        /* Remove box shadow */
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }}
        .logo:hover {{
            transform: translateY(-3px);
            box-shadow: none;        /* No shadow on hover */
            background: transparent; /* Keep transparent background on hover */
        }}
        .logo img {{
            width: 220px;
            height: auto;
            border-radius: 8px;
            background: transparent;   /* Ensure image background is transparent */
            box-shadow: none;          /* Remove any shading */
        }}
        .sidebar-section {{
            background: rgba(255, 255, 255, 0.95);
            margin: 15px;
            padding: 20px;
            border-radius: 16px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
            border: 1px solid rgba(255, 255, 255, 0.3);
            backdrop-filter: blur(10px);
        }}
        .sidebar-title {{
            font-size: 16px;
            font-weight: 700;
            color: #2d3748;
            margin-bottom: 15px;
            text-align: center;
            background: linear-gradient(135deg, #667eea, #764ba2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            padding-bottom: 8px;
            border-bottom: 2px solid #e2e8f0;
        }}
        .sidebar-widget {{
            background: linear-gradient(135deg, #f8fafc, #f1f5f9);
            padding: 15px;
            border-radius: 12px;
            margin: 10px 0;
            border: 1px solid #e2e8f0;
            transition: all 0.3s ease;
        }}
        .sidebar-widget:hover {{
            transform: translateX(5px);
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.2);
            border-color: #667eea;
        }}
        [data-testid="stSidebar"] {{
            background: linear-gradient(180deg, #667eea 0%, #764ba2 100%);
        }}
        [data-testid="stSidebar"] * {{
            color: white;
        }}
        /* Fix: Override text input and textarea text color for sidebar */
        [data-testid="stSidebar"] input,
        [data-testid="stSidebar"] textarea {{
            color: #222 !important;
            background: rgba(255, 255, 255, 0.92);
        }}
        [data-testid="stSidebar"] input::placeholder,
        [data-testid="stSidebar"] textarea::placeholder {{
            color: #667eea99 !important;
        }}
        /* Style the Streamlit radio widget on the sidebar to match sidebar theme */
        [data-testid="stSidebar"] [data-baseweb="radio"] label {{
            color: black !important;
            background: transparent;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] .css-1aehpvj,  /* radio background */
        [data-testid="stSidebar"] [data-baseweb="radio"] .css-1e3y7mc {{
            background: black !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] span {{
            color: white !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] svg {{
            fill: #e7e5f8 !important;
            stroke: #fff !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] input[type="radio"]:checked + div > div > svg {{
            fill: #fff !important;
            stroke: #764ba2 !important;
        }}
        </style>

        """,
        unsafe_allow_html=True,
    )

        
    # <div class="logo-container">
    #     <div class="logo">
    #         <img src="data:image/png;base64,{logo_base64}" alt="Logo">
    #     </div>
    # </div>
    main()
    

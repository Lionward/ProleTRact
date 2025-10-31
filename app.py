import streamlit as st
import streamlit.components.v1 as components
from proletract.modules.io.data_handling import VCFHandler, CohortHandler
from proletract.modules.viz.visualization import Visualization
import pandas as pd 
import base64
import os
from pathlib import Path



def get_pathogenic_TRs():
    # Resolve to package data file
    root = Path(__file__).resolve().parent
    data_path = root / "proletract" / "data" / "pathogenic_TRs.bed"
    pathogenic_trs = pd.read_csv(str(data_path), sep="\t", header=None)
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
        ("individual sample 👤", "Cohort 👥👥", "comparison 🔄"), 
        key="analysis_mode_radio",
        label_visibility='visible',
        help="Choose the analysis workflow.",
    )
    st.markdown("""
        <style>
            /* Radio button container with white/light background to contrast with purple sidebar */
            [data-testid="stSidebar"] div[role="radiogroup"] {
                background: rgba(255, 255, 255, 0.95) !important;
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1) !important;
                border: 1px solid rgba(255, 255, 255, 0.3) !important;
                border-radius: 12px !important;
                padding: 12px 10px !important;
                margin: 10px 0 !important;
                backdrop-filter: blur(10px);
            }
            /* Make internal radio controls and labels transparent/native and styled */
            [data-testid="stSidebar"] div[data-baseweb="radio"] {
                background: transparent !important;
                box-shadow: none !important;
            }
            [data-testid="stSidebar"] div[data-baseweb="radio"] label {
                font-weight: 600 !important;
                color: #1f2937 !important;
                background: transparent !important;
                padding: 8px 12px !important;
                border-radius: 8px !important;
                transition: all 0.2s ease !important;
            }
            /* Force text color to be dark */
            [data-testid="stSidebar"] div[data-baseweb="radio"] label,
            [data-testid="stSidebar"] div[data-baseweb="radio"] label span,
            [data-testid="stSidebar"] div[data-baseweb="radio"] label * {
                color: #1f2937 !important;
            }
            [data-testid="stSidebar"] div[data-baseweb="radio"] label:hover {
                background: rgba(102, 126, 234, 0.1) !important;
            }
            [data-testid="stSidebar"] div[data-baseweb="radio"] input[type="radio"]:checked + div label {
                background: rgba(102, 126, 234, 0.15) !important;
            }
            [data-testid="stSidebar"] div[data-baseweb="radio"] input[type="radio"]:checked + div label,
            [data-testid="stSidebar"] div[data-baseweb="radio"] input[type="radio"]:checked + div label span,
            [data-testid="stSidebar"] div[data-baseweb="radio"] input[type="radio"]:checked + div label * {
                color: #4f46e5 !important;
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


    if analysis_mode == "individual sample 👤":

        vcf_handler.handle_individual_sample()
        st.session_state.analysis_mode = "individual sample"
        
        if 'records_map' in st.session_state:
            
            visualization.visulize_region() 
        else:
            st.stop()


    elif analysis_mode == "Cohort 👥👥":
        # Creative toggle-style selector for cohort mode
        # Initialize cohort_mode if not exists
        if 'cohort_mode' not in st.session_state:
            st.session_state.cohort_mode = "reads"
        

        # Creative toggle-style selector using Streamlit buttons styled like cards
        
        reads_selected = st.session_state.cohort_mode == "reads"
        assembly_selected = st.session_state.cohort_mode == "assembly"
        
        # --- Custom nicer style for toggle buttons ---
        st.sidebar.markdown("""
            <style>
            .cohort-toggle-btn {
                display: block;
                width: 100%;
                background: linear-gradient(93deg, #667eea 0%, #764ba2 100%);
                color: white !important;
                border: none;
                border-radius: 14px;
                padding: 18px 16px 15px 16px;
                margin-bottom: 10px;
                font-size: 17px !important;
                font-weight: 800 !important;
                box-shadow: 0 2px 14px 0 rgba(102,126,234,0.14);
                letter-spacing: 0.01em;
                text-align: left;
                outline: none;
                transition: all .2s;
                position: relative;
                cursor: pointer;
                line-height: 1.22;
            }
            .cohort-toggle-btn.selected {
                background: linear-gradient(99deg, #764ba2 6%, #667eea 94%);
                color: #fff;
                box-shadow: 0 3px 38px 0 rgba(118, 75, 162, 0.18);
                border: 2.2px solid #fff5;
                opacity: 1.0;
                transform: scale(1.034);
            }
            .cohort-toggle-btn:hover {
                filter: brightness(1.08);
                transform: translateY(-2.5px) scale(1.04);
                background: linear-gradient(99deg, #7b88ee 10%, #a38fcf 100%);
            }
            .cohort-toggle-btn .subtitle {
                font-size: 13px !important;
                font-weight: 500;
                opacity: 0.85;
                display: block;
                margin-top: 0.3em;
                letter-spacing: 0;
                color: #f3f1ff;
            }
            </style>
        """, unsafe_allow_html=True)
        

        # Ensure cohort_mode is initialized
        if st.session_state.get('cohort_mode') not in ("reads", "assembly"):
            st.session_state['cohort_mode'] = "reads"
        
        # Get current mode
        current_mode = st.session_state.get('cohort_mode', 'reads')
        reads_selected = current_mode == "reads"
        assembly_selected = current_mode == "assembly"
        
        # Create buttons side by side
        col1, col2 = st.sidebar.columns(2)
        
        with col1:
            # Add indicator to reads button if active
            reads_label = "☰ Reads-based VCF" + (" ✓" if reads_selected else "")
            reads_btn_clicked = st.button(
                reads_label, 
                key="reads_btn", 
                help="Per individual TR genotyping",
                use_container_width=True,
                type="primary" if reads_selected else "secondary"
            )
            if reads_btn_clicked:
                st.session_state['cohort_mode'] = "reads"
                st.rerun()
        
        with col2:
            # Add indicator to assembly button if active
            assembly_label = "━━ Assembly VCF" + (" ✓" if assembly_selected else "")
            assembly_btn_clicked = st.button(
                assembly_label, 
                key="assembly_btn", 
                help="Haplotype-resolved",
                use_container_width=True,
                type="primary" if assembly_selected else "secondary"
            )
            if assembly_btn_clicked:
                st.session_state['cohort_mode'] = "assembly"
                st.rerun()

        cohort_handler.handle_cohort()
        st.session_state.analysis_mode = "Cohort"
        visualization.visulize_cohort()
        # if the path doesn't end with a slash add it
    elif analysis_mode == "comparison 🔄":
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
    with st.sidebar.expander("✨ About ProleTRact", expanded=False):
        st.info(
            "ProleTRact is a tool for exploring, visualizing, and comparing tandem repeat regions using VCF files from TandemTwister output. Switch modes below to begin.",
            icon="💡"
        )
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
        /* Exception: Radio buttons in white container should have dark text */
        [data-testid="stSidebar"] div[role="radiogroup"] * {{
            color: #1f2937 !important;
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
        /* Style the Streamlit radio widget on the sidebar - ensure text is visible on white background */
        [data-testid="stSidebar"] [data-baseweb="radio"] label {{
            color: #1f2937 !important;
            background: transparent;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] .css-1aehpvj,  /* radio background */
        [data-testid="stSidebar"] [data-baseweb="radio"] .css-1e3y7mc {{
            background: #667eea !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] span {{
            color: #1f2937 !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] label span,
        [data-testid="stSidebar"] [data-baseweb="radio"] label * {{
            color: #1f2937 !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] input[type="radio"]:checked + div label,
        [data-testid="stSidebar"] [data-baseweb="radio"] input[type="radio"]:checked + div label span,
        [data-testid="stSidebar"] [data-baseweb="radio"] input[type="radio"]:checked + div label * {{
            color: #4f46e5 !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] svg {{
            fill: #667eea !important;
            stroke: #fff !important;
        }}
        [data-testid="stSidebar"] [data-baseweb="radio"] input[type="radio"]:checked + div > div > svg {{
            fill: #667eea !important;
            stroke: #fff !important;
        }}
        </style>

        """,
        unsafe_allow_html=True,
    )

    main()
    

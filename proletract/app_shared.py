import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path
import pandas as pd
from proletract.modules.io.data_handling import VCFHandler, CohortHandler
from proletract.modules.viz.visualization import Visualization


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


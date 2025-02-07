import streamlit as st
from data_handling import VCFHandler, CohortHandler
from visualization import Visualization
import pandas as pd 
import base64


def get_image_base64(path):

    with open(path, "rb") as f:
        encoded_string = base64.b64encode(f.read()).decode()
    return encoded_string

def get_pathogenic_TRs():
    pathogenic_trs = pd.read_csv("/confidential/tGenVar/Lion/TandemTwist/tr_regions/pathogeic_TRs.bed", sep="\t", header=None)
    pathogenic_trs.columns = ["chrom",'start','end','motif','pathogenic_min','inheritance','disease','gene']
    st.session_state.pathogenic_TRs = pathogenic_trs
    # make column region
    pathogenic_trs['region'] = pathogenic_trs['chrom'] + ":" + pathogenic_trs['start'].astype(str) + "-" + pathogenic_trs['end'].astype(str)


def main():

    if 'vcf_handler' not in st.session_state:
        st.session_state.vcf_handler = VCFHandler()
    if 'cohort_handler' not in st.session_state:
        st.session_state.cohort_handler = CohortHandler()
    if 'visualization' not in st.session_state:
        st.session_state.visualization = Visualization()
    get_pathogenic_TRs()
    
    vcf_handler = st.session_state.vcf_handler
    cohort_handler = st.session_state.cohort_handler
    visualization = st.session_state.visualization

    #st.session_state.pathogenic_TRs = PATHOGENIC_TRS

    st.sidebar.markdown("### Tandem Repeat Visualization")
    analysis_mode = st.sidebar.radio("Select the type of analysis", ("individual sample", "Cohort"))

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




if __name__ == "__main__":
    st.set_page_config(layout="wide")
    # Path to your local logo
    logo_path = "tandem_twister_vis_logo.png"
    logo_base64 = get_image_base64(logo_path)
    # Custom HTML/CSS to align the logo to the left with a transparent background
    st.sidebar.markdown(
        f"""
        <style>
        .logo {{
            display: flex;
            justify-content: center;
            align-items: center;
            margin-bottom: 20px;
            background-color: transparent;  /* Ensure background is transparent */
        }}
        .logo img {{
            width: 250px;  /* Adjust width as needed */
        }}
        </style>
        <div class="logo">
            <img src="data:image/png;base64,{logo_base64}" alt="Logo">
        </div>
        """,
        unsafe_allow_html=True,
    )

    
    
    main()
    

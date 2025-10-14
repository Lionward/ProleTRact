import pysam
import os
import re
import streamlit as st

class VCFHandler:
    def __init__(self):
        self.vcf_file_path = None
        self.records = st.session_state.get('records', None)
        self.records_map = st.session_state.get('records_map', None)

    def parse_vcf(self, vcf_file):

        vcf = pysam.VariantFile(vcf_file)
        records_ids = {}
        records_map = {}
        idx = 0
        for rec in vcf.fetch():
            records_ids[rec.id] = f"{rec.chrom}:{rec.pos}-{rec.stop}"
            records_map[idx] = rec.id
            idx += 1
        return records_ids, records_map


    def load_vcf(self,vcf_file):
        return pysam.VariantFile(vcf_file)
    

    def handle_individual_sample(self):

        if 'vcf_file_path' not in st.session_state:
            # initialize the vcf file path
            st.session_state.vcf_file_path = None
        vcf_path = st.sidebar.text_input("Enter the path of your VCF file", key="vcf_file_path_input", help="Enter the path of your VCF file, the file should be zipped and indexed with tabix")
        _, _, middle, _ = st.sidebar.columns([1, 0.3, 2, 1])

        if middle.button("Upload VCF File"):
            if vcf_path:
                st.session_state.vcf_file_path = vcf_path
                st.session_state.pop('records', None)
                st.session_state.pop('records_map', None)

            else:
                st.info("Please enter the path to the VCF file")
            if 'records' not in st.session_state:
                if st.session_state.vcf_file_path:
                    st.session_state.records, st.session_state.records_map = self.parse_vcf( st.session_state.vcf_file_path)
                    st.session_state.hgsvc_path = "/confidential/Structural_Variation/output_maryam/pipeline/input_folder/asm_samples/assembly_results/"
                    # check if the path exists
                    if os.path.exists(st.session_state.hgsvc_path):
                        st.session_state.file_paths = [f for f in os.listdir(st.session_state.hgsvc_path) if f.endswith('h1.vcf.gz') or f.endswith('h2.vcf.gz')]
                        st.session_state.files = [self.load_vcf(st.session_state.hgsvc_path + f) for f in st.session_state.file_paths]
                    else:
                        st.session_state.files = None
                        st.session_state.file_paths = None
                else:
                    st.error("VCF file path is not set.")

class CohortHandler(VCFHandler):
    def __init__(self):
            # self.cohort_path = st.session_state.get('cohort_path', None)
        pass
    def handle_cohort(self):
         
        if 'path_to_cohort' not in st.session_state:
            st.session_state.path_to_cohort = None
        cohort_path = st.sidebar.text_input("Enter the path to the cohort results", key="cohort_path_input", help="Enter the path to the cohort results, the files should be zipped and indexed with tabix")

        if cohort_path is None:
            st.stop()
        if not cohort_path.endswith('/'):
            cohort_path += '/'
            

        if st.sidebar.button("Load Cohort"):
            if cohort_path:
                st.session_state.path_to_cohort = cohort_path
            st.session_state.cohort_file_paths = [f for f in os.listdir(st.session_state.path_to_cohort) if f.endswith('.vcf.gz')]
            st.session_state.cohort_files = [self.load_vcf(st.session_state.path_to_cohort + f) for f in st.session_state.cohort_file_paths]
            
            st.session_state.cohorts_records_map = self.get_records_info(st.session_state.path_to_cohort + st.session_state.cohort_file_paths[0])

    def get_records_info(self, vcf_file):
        vcf = pysam.VariantFile(vcf_file)
        cohorts_map = {}
        idx = 0
        for rec in vcf:
            cohorts_map[idx] = rec.id
            idx += 1
        return cohorts_map
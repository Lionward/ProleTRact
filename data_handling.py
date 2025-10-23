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

    def parse_vcf_TRGT(self, vcf_file):
        vcf = pysam.VariantFile(vcf_file)
        records_map = {}
        records_ids = {}
        idx = 0
        for rec in vcf.fetch():
            trid = rec.info['TRID']
            chrom, start, end = trid.split('_')
            records_ids[f"{chrom}:{start}-{end}"] = f"{rec.chrom}:{rec.pos}-{rec.stop}"
            records_map[idx] = f"{rec.chrom}:{rec.pos}-{rec.stop}"
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
        with st.spinner("Wait for it..."):
            button_clicked = middle.button(
                "Upload VCF File",
                key="upload_vcf_btn",
                help=None,
                type="secondary",
                use_container_width=False,
                kwargs={
                    "style": "font-size: 12px !important; padding: 4px 16px !important;"
                }
            )
            if button_clicked:
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
            
        
    def handle_comparison_samples(self):
        if 'vcf_file_path_input1' not in st.session_state:
            st.session_state.vcf_file_path_input1 = None
        if 'vcf_file_path_input2' not in st.session_state:
            st.session_state.vcf_file_path_input2 = None
        if 'assembly_vcf_input_h1' not in st.session_state:
            st.session_state.assembly_vcf_input_h1 = None
        if 'assembly_vcf_input_h2' not in st.session_state:
            st.session_state.assembly_vcf_input_h2 = None
        # Set default values for the comparison input fields
        default_trgt = "/confidential/home01/Calraei/tandemrepeats/trgt/results/HG002_trgt_sorted.vcf.gz"
        default_tandemtwister = "/confidential/home01/Calraei/tandemrepeats/tandemtwister/results/HG002_perfect_CCS_tandemtwister.vcf.gz"
        default_assembly_h1 = "/confidential/home01/Calraei/tandemrepeats/tandemtwister/assembly_results/NA24385_h1.vcf.gz"
        default_assembly_h2 = "/confidential/home01/Calraei/tandemrepeats/tandemtwister/assembly_results/NA24385_h2.vcf.gz"

        # initialize the input fields in session state if not yet present
        if 'vcf_file_trgt_input' not in st.session_state:
            st.session_state.vcf_file_trgt_input = default_trgt
        if 'vcf_file_tandemtwister_input' not in st.session_state:
            st.session_state.vcf_file_tandemtwister_input = default_tandemtwister
        if 'assembly_vcf_h1_input' not in st.session_state:
            st.session_state.assembly_vcf_h1_input = default_assembly_h1
        if 'assembly_vcf_h2_input' not in st.session_state:
            st.session_state.assembly_vcf_h2_input = default_assembly_h2

        vcf_file_trgt = st.sidebar.text_input(
            "Enter the path of your TRGT VCF file",
            key="vcf_file_trgt_input",
            help="Enter the path of your TRGT VCF file, the file should be zipped and indexed with tabix",
            value=st.session_state.vcf_file_trgt_input
        )
        vcf_file_tandemtwister = st.sidebar.text_input(
            "Enter the path of your TandemTwister VCF file",
            key="vcf_file_tandemtwister_input",
            help="Enter the path of your TandemTwister VCF file, the file should be zipped and indexed with tabix",
            value=st.session_state.vcf_file_tandemtwister_input
        )
        assembly_vcf_h1 = st.sidebar.text_input(
            "Enter the path of your assembly VCF file",
            key="assembly_vcf_h1_input",
            help="Enter the path of your assembly haplotype 1 VCF file, the file should be zipped and indexed with tabix",
            value=st.session_state.assembly_vcf_h1_input
        )
        assembly_vcf_h2 = st.sidebar.text_input(
            "Enter the path of your assembly VCF file",
            key="assembly_vcf_h2_input",
            help="Enter the path of your assembly haplotype 2 VCF file, the file should be zipped and indexed with tabix",
            value=st.session_state.assembly_vcf_h2_input
        )
        _, _, middle, _ = st.sidebar.columns([1, 0.3, 2, 1])

        if middle.button("Load Comparison Samples"):
            if vcf_file_trgt and vcf_file_tandemtwister and assembly_vcf_h1 and assembly_vcf_h2:
                st.session_state.vcf_file_trgt = vcf_file_trgt
                st.session_state.vcf_file_tandemtwister = vcf_file_tandemtwister
                st.session_state.assembly_vcf_h1 = assembly_vcf_h1
                st.session_state.assembly_vcf_h2 = assembly_vcf_h2
                st.session_state.pop('vcf_file_trgt_records_map', None)
                st.session_state.pop('vcf_file_tandemtwister_records_map', None)
                st.session_state.pop('assembly_vcf_records_map', None)
                st.session_state.pop('vcf_file_trgt_records', None)
                st.session_state.pop('vcf_file_tandemtwister_records', None)
                st.session_state.pop('assembly_vcf_records_h1', None)
                st.session_state.pop('assembly_vcf_records_h2', None)

            else:  
                st.info("Please enter the path to the VCF files")
                st.stop()
            if 'vcf_file_trgt_records' not in st.session_state:
                try:
                    st.session_state.vcf_file_trgt_records, st.session_state.vcf_file_trgt_records_map = self.parse_vcf_TRGT(st.session_state.vcf_file_trgt)
                except:
                    info_container.warning("Failed to parse the TRGT VCF file")
            if 'vcf_file_tandemtwister_records' not in st.session_state:
                info_container = st.container()
                try:
                    st.session_state.vcf_file_tandemtwister_records, st.session_state.vcf_file_tandemtwister_records_map = self.parse_vcf(st.session_state.vcf_file_tandemtwister)
                except:
                    info_container.warning("Failed to parse the TandemTwister VCF file")
            if 'assembly_vcf_records_h1' not in st.session_state:
                info_container = st.container()
                try:
                    st.session_state.assembly_vcf_records_h1, st.session_state.assembly_vcf_records_map_h1 = self.parse_vcf(st.session_state.assembly_vcf_h1)
                except:
                    info_container.warning("Failed to parse the assembly VCF file for haplotype 1")
            if 'assembly_vcf_records_h2' not in st.session_state:
                info_container = st.container()
                try:
                    st.session_state.assembly_vcf_records_h2, st.session_state.assembly_vcf_records_map_h2 = self.parse_vcf(st.session_state.assembly_vcf_h2)
                except:
                    info_container.warning("Failed to parse the assembly VCF file for haplotype 2")
            st.session_state.all_files_parsed = True
        elif st.session_state.all_files_parsed and st.session_state.vcf_file_trgt_records and st.session_state.vcf_file_tandemtwister_records and st.session_state.assembly_vcf_records_h1 and st.session_state.assembly_vcf_records_h2:
            pass
        else:
            st.stop()

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


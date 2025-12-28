import pysam
import os
import re
import streamlit as st

# Try to import Cython-optimized functions, fallback to pure Python
try:
    from .fast_parsing import parse_vcf_record_ids_fast as _parse_vcf_record_ids_fast
    _USE_CYTHON = True
    print("🚀 Cython optimizations: ENABLED (fast_parsing module loaded)")
except (ImportError, ModuleNotFoundError):
    _USE_CYTHON = False
    _parse_vcf_record_ids_fast = None
    print("⚠️  Cython optimizations: DISABLED (using pure Python fallback)")

class VCFHandler:
    def __init__(self):
        self.vcf_file_path = None
        self.records = st.session_state.get('records', None)
        self.records_map = st.session_state.get('records_map', None)
        # Initialize region_genotypes if not exists (backward compatibility)
        if 'region_genotypes' not in st.session_state:
            st.session_state.region_genotypes = {}

    def parse_vcf(self, vcf_file):
        """
        Parse VCF file and extract region IDs, mapping, and genotypes.
        Returns (records_ids, records_map, region_genotypes)
        where region_genotypes maps region string -> genotype string
        """
        vcf = pysam.VariantFile(vcf_file)
        
        # Try Cython version first (now includes genotype extraction)
        if _USE_CYTHON and _parse_vcf_record_ids_fast is not None:
            try:
                records_ids, records_map, region_genotypes = _parse_vcf_record_ids_fast(vcf)
                return records_ids, records_map, region_genotypes
            except Exception:
                # Fallback to pure Python if Cython fails
                pass
        
        # Python version with genotype extraction
        records_ids = {}  # dict for mapping region_str -> record_id
        records_map = {} # dict for mapping idx -> record_id
        region_genotypes = {}  # dict for mapping region_str -> genotype
        idx = 0
        for rec in vcf.fetch():
            region_str = f"{rec.chrom}:{rec.pos}-{rec.stop}"
            records_ids[region_str] = rec.id  #
            records_map[idx] = rec.id
            
            # Extract genotype
            try:
                gt = rec.samples[0]['GT']
                if gt is not None:
                    gt_str = '/'.join([str(i) for i in gt]) if isinstance(gt, (tuple, list)) else str(gt)
                    region_genotypes[region_str] = gt_str
                else:
                    region_genotypes[region_str] = './.'
            except (KeyError, IndexError, AttributeError):
                region_genotypes[region_str] = './.'
            
            idx += 1
        return records_ids, records_map, region_genotypes



    def load_vcf(self,vcf_file):
        return pysam.VariantFile(vcf_file)
    

    def handle_individual_sample(self):

        if 'vcf_file_path' not in st.session_state:
            # initialize the vcf file path
            st.session_state.vcf_file_path = None
        with st.sidebar.expander("📂 input data", expanded=True):
            vcf_path = st.text_input(
                "Enter the path of your VCF file",
                key="vcf_file_path_input",
                help="Enter the path of your VCF file, the file should be zipped and indexed with tabix",
            )

            public_vcf_folder = st.text_input(
                "Enter the path of the public VCF folder",
                key="public_vcf_folder_input",
                help="Enter the path of the public VCF folder, the folder should contain the VCF files",
            )
            if not public_vcf_folder.endswith('/'):
                public_vcf_folder += '/'
                
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
                    st.session_state.pop('read_support', None)
                    st.session_state.pop('read_support_source', None)
                else:
                    st.info("Please enter the path to the VCF file")

                if 'records' not in st.session_state:
                    if st.session_state.vcf_file_path:
                        # Use cached parsing with progress indicator
                        with st.spinner("Parsing VCF file... This may take a moment."):
                            st.session_state.records, st.session_state.records_map, st.session_state.region_genotypes = self.parse_vcf(st.session_state.vcf_file_path)
                        st.session_state.hgsvc_path = public_vcf_folder 
                        # check if the path exists
                        if os.path.exists(st.session_state.hgsvc_path):
                            with st.spinner("Loading public VCF files..."):
                                st.session_state.file_paths = [f for f in os.listdir(st.session_state.hgsvc_path) if f.endswith('h1.vcf.gz') or f.endswith('h2.vcf.gz')]
                                st.session_state.files = [self.load_vcf(st.session_state.hgsvc_path + f) for f in st.session_state.file_paths]
                        else:
                            st.session_state.files = None
                            st.session_state.file_paths = None
                    else:
                        st.error("VCF file path is not set.")
            

        
class CohortHandler(VCFHandler):
    def __init__(self):
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
            with st.spinner("Loading cohort files..."):
                st.session_state.cohort_file_paths = [f for f in os.listdir(st.session_state.path_to_cohort) if f.endswith('.vcf.gz')]
                st.session_state.cohort_files = [self.load_vcf(st.session_state.path_to_cohort + f) for f in st.session_state.cohort_file_paths]
            
            with st.spinner("Parsing cohort records..."):
                st.session_state.cohorts_records_map, st.session_state.cohort_region_genotypes = self.get_records_info(st.session_state.path_to_cohort + st.session_state.cohort_file_paths[0])
    
    def get_records_info(self, vcf_file):
        """
        Get records mapping and genotype information for cohort.
        Returns (cohorts_map, region_genotypes) where region_genotypes maps region_str -> genotype
        Note: Uses record.id as key to match cohorts_map values
        """
        vcf = pysam.VariantFile(vcf_file)
        cohorts_map = {}
        region_genotypes = {}
        idx = 0
        for rec in vcf:
            cohorts_map[idx] = rec.id
            
            # Extract genotype from first sample - use rec.id as key to match cohorts_map
            try:
                gt = rec.samples[0]['GT']
                if gt is not None:
                    gt_str = '/'.join([str(i) for i in gt]) if isinstance(gt, (tuple, list)) else str(gt)
                    region_genotypes[rec.id] = gt_str
                else:
                    region_genotypes[rec.id] = './.'
            except (KeyError, IndexError, AttributeError):
                region_genotypes[rec.id] = './.'
            
            idx += 1
        return cohorts_map, region_genotypes


"""
Display functions for visualization.
"""
import streamlit as st
import pandas as pd
import re
from proletract.modules.viz import vis_helper as _vh
motif_legend_html = _vh.motif_legend_html
plot_motif_bar = _vh.plot_motif_bar
display_motifs_as_bars = _vh.display_motifs_as_bars
display_genotype_card = _vh.display_genotype_card
display_dynamic_sequence_with_highlighted_motifs = _vh.display_dynamic_sequence_with_highlighted_motifs
create_genotype_comparison_matrix = _vh.create_genotype_comparison_matrix
from proletract.modules.viz import parsers
from proletract.modules.viz import plots
from proletract.modules.viz import utils


def render_region_display(markdown_placeholder, region):
        # Parse region to create links to various genomic databases
        chrom = None
        pos_range = None
        chrom_no_chr = None
        start = None
        end = None
        
        try:
            if ':' in region:
                chrom, pos_range = region.split(':', 1)
                chrom_no_chr = chrom.replace('chr', '').replace('CHR', '')
                
                # Parse start and end positions
                if '-' in pos_range:
                    start, end = pos_range.split('-', 1)
                    start = start.strip()
                    end = end.strip()
                
                # Ensure chromosome format is correct (add 'chr' if needed)
                if not chrom.startswith('chr'):
                    chrom = f'chr{chrom}'
        except Exception:
            pass
        
        # Generate URLs for various databases
        urls = {}
        if chrom and pos_range:
            # UCSC Genome Browser
            urls['UCSC'] = {
                'url': f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}:{pos_range}",
                'icon': '🌐',
                'color': 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)'
            }
            
            # Ensembl (for human, use GRCh38)
            if chrom_no_chr:
                urls['Ensembl'] = {
                    'url': f"https://www.ensembl.org/Homo_sapiens/Location/View?r={chrom_no_chr}:{pos_range}",
                    'icon': '🧬',
                    'color': 'linear-gradient(135deg, #FF6B6B 0%, #EE5A6F 100%)'
                }
            
            # NCBI Genome Data Viewer
            if chrom_no_chr:
                urls['NCBI'] = {
                    'url': f"https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000001405.40&chr={chrom_no_chr}&from={start}&to={end}" if start and end else f"https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000001405.40&chr={chrom_no_chr}",
                    'icon': '📊',
                    'color': 'linear-gradient(135deg, #96CEB4 0%, #6C9A8B 100%)'
                }
            
            # gnomAD (for variant frequency data)
            if chrom_no_chr and start and end:
                urls['gnomAD'] = {
                    'url': f"https://gnomad.broadinstitute.org/region/{chrom_no_chr}-{start}-{end}",
                    'icon': '📈',
                    'color': 'linear-gradient(135deg, #FFE66D 0%, #FF6B6B 100%)'
                }
            
            # DECIPHER (for pathogenic variants and phenotype data)
            if chrom_no_chr:
                urls['DECIPHER'] = {
                    'url': f"https://www.deciphergenomics.org/browser#q/grch37:{chrom_no_chr}:{pos_range}",
                    'icon': '🔍',
                    'color': 'linear-gradient(135deg, #A8E6CF 0%, #7FCDCD 100%)'
                }
            
            # TRExplorer - Tandem Repeat Explorer catalog
            if chrom and pos_range:
                # TRExplorer uses searchQuery parameter with chr:start-end format
                # Format: https://trexplorer.broadinstitute.org/#searchQuery=chr1:94418422-94418444
                urls['TRExplorer'] = {
                    'url': f"https://trexplorer.broadinstitute.org/#sc=isPathogenic&sd=DESC&showRs=1&searchQuery={chrom}:{pos_range}&showColumns=0i1i2i3i4i7i21i17",
                    'icon': '🔬',
                    'color': 'linear-gradient(135deg, #10B981 0%, #059669 100%)'
                }
        
        html = f"""
            <style>
                .region-badge {{
                    display: inline-block;
                    background: #ECEAFB;
                    padding: 6px 12px;
                    border-radius: 50px;
                    font-size: 8px;
                    font-weight: 700;
                    margin-bottom: 20px;
                    backdrop-filter: blur(10px);
                    border: 1px solid #C9BEEF;
                    letter-spacing: 1px;
                }}
                .external-link {{
                    display: inline-flex;
                    align-items: center;
                    gap: 10px;
                    color: white !important;
                    text-decoration: none;
                    padding: 18px 28px;
                    border-radius: 16px;
                    font-size: 20px;
                    font-weight: 700;
                    transition: all 0.3s ease;
                    box-shadow: 0 3px 10px rgba(0, 0, 0, 0.22);
                }}
                .external-link:hover {{
                    transform: translateY(-4px) scale(1.04);
                    box-shadow: 0 6px 16px rgba(0, 0, 0, 0.32);
                }}
                .external-link-icon {{
                    font-size: 28px;
                }}
                .links-container {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 18px;
                    justify-content: center;
                    margin-top: 23px;
                }}
            </style>
            
            <div style="display: flex; flex-direction: column; justify-content: center; align-items: center; min-height: 80px;">
                <div style="
                    display: flex; 
                    align-items: center; 
                    background: #ECEAFB; 
                    padding: 13px 24px; 
                    border-radius: 24px; 
                    box-shadow: 0 4px 16px 2px rgba(118, 75, 162, 0.12); 
                    font-size: 2rem; 
                    font-weight: 700;
                    border: 1px solid #C9BEEF;
                    margin-bottom: 15px;
                    ">
                    <span style="width: 26px; height: 26px; background: #764ba2; border-radius: 50%; margin-right: 28px; box-shadow: 0 0 12px 2px #A184D6;"></span>
                    Region: {region}
                </div>
                <div class="links-container">
                    {''.join([f'<a href="{info["url"]}" target="_blank" class="external-link" style="background: {info["color"]};"><span class="external-link-icon">{info["icon"]}</span> {name}</a>' for name, info in urls.items()])}
                </div>
            </div>
        """
        markdown_placeholder.html(html)


def display_motifs_with_bars(record, containter, motif_colors, CN1_col, CN2_col, show_comparison, hgsvc_records):
        motif_names, motif_count_h1, motif_count_h2 = parsers.parse_motif_in_region(record, parsers.count_motifs)
        if motif_count_h1 == None and motif_count_h2 == None:
            st.info("No motifs found in the region")
            return
        st.markdown("""
            <style>
                .stTabs [data-baseweb="tab-list"] {
                    gap: 6px;
                    background: #f8fafc;
                    padding: 4px 4px;
                    border-radius: 16px;
                    margin-bottom: 18px;
                    min-height: 22px;
                }
                .stTabs [data-baseweb="tab"] {
                    height: 28px;
                    min-height: 28px;
                    min-width: 38px;
                    white-space: pre-wrap;
                    background: white;
                    border-radius: 10px;
                    border: 2px solid #e2e8f0;
                    color: #64748b;
                    font-weight: 700;
                    font-size: 1.13rem;
                    padding: 0 13px;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                    line-height: 1.33;
                }
                .stTabs [data-baseweb="tab"]:hover {
                    background: #f1f5f9;
                    border-color: #cbd5e1;
                    color: #475569;
                    transform: translateY(-1.8px) scale(1.04);
                    box-shadow: 0 2px 12px rgba(0,0,0,0.13);
                }
                .stTabs [aria-selected="true"] {
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                    color: white !important;
                    border-color: #667eea !important;
                    box-shadow: 0 5px 19px rgba(102, 126, 234, 0.23) !important;
                    transform: translateY(-2px) scale(1.045);
                    font-size: 1.23rem !important;
                }
                .stTabs [data-baseweb="tab-highlight"] {
                    background-color: transparent !important;
                }
            </style>
        """, unsafe_allow_html=True)

        

        with containter:
            tab1, tab2, tab3 = st.tabs([
                "🧬 **Alleles**", 
                "🔄 **Alleles vs Ref**", 
                "🌍 **Alleles vs Pop**"
            ])
            st.markdown(
                "<style>.tab-content {font-size: 20px;}</style>",
                unsafe_allow_html=True,
            )


            with tab2:
                motif_legend_html(record['motif_ids_ref'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("### Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                display_motifs_as_bars("Ref", motif_colors, record['motif_ids_ref'], record['spans'][0], record['ref_allele'], motif_names, None)
                display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names, record['supported_reads_h1'])
                if record['alt_allele2'] != '':
                    display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names, record['supported_reads_h2'])
            with tab1:
                motif_legend_html(record['motif_ids_h1'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("###  Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names, record['supported_reads_h1'])
                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        with st.expander("Show motif bar plot for Allel1", expanded=False):
                            plot_motif_bar(motif_count_h1, motif_names, motif_colors)
                
                if record['alt_allele2'] != '':
                    display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names, record['supported_reads_h2'])
                    plot_container_h2 = st.empty()

                    with plot_container_h2:
                        if show_comparison == False:
                            with st.expander("Show motif bar plot for Allel2", expanded=False):
                                plot_motif_bar(motif_count_h2, motif_names, motif_colors)

            with tab3:
                if hgsvc_records:
                    # Add genotype comparison for population data
                    # For current sample: compute diploid genotype from haplotypes (similar to cohort mode)
                    current_sample_gt = record['gt']
                    # Extract gt_h1 and gt_h2 from the GT field if available
                    # The GT field format is "0/1" or "0/0", so we need to split it
                    gt_parts = current_sample_gt.split('/') if '/' in str(current_sample_gt) else [current_sample_gt, current_sample_gt]
                    gt_h1 = gt_parts[0] if len(gt_parts) > 0 else "0"
                    gt_h2 = gt_parts[1] if len(gt_parts) > 1 else gt_parts[0]
                    
                    # Compute diploid genotype using the same logic as cohort mode
                    diploid_gt = parsers.compute_diploid_genotype_assembly(
                        gt_h1, gt_h2, 
                        record.get('motif_ids_h1', []), 
                        record.get('motif_ids_h2', [])
                    )
                    
                    population_genotypes = {"Current Sample": diploid_gt}
                    
                    # Handle population records - group by base name if they have _h1/_h2 suffixes
                    sample_groups = {}
                    for sample_name, sample_record in hgsvc_records.items():
                        if sample_name.endswith('_h1'):
                            base_name = sample_name[:-3]
                            if base_name not in sample_groups:
                                sample_groups[base_name] = {}
                            sample_groups[base_name]['h1'] = sample_record
                        elif sample_name.endswith('_h2'):
                            base_name = sample_name[:-3]
                            if base_name not in sample_groups:
                                sample_groups[base_name] = {}
                            sample_groups[base_name]['h2'] = sample_record
                        else:
                            # Sample doesn't have _h1/_h2 suffix, use as is
                            if 'gt' in sample_record:
                                population_genotypes[sample_name] = sample_record['gt']
                    
                    # Compute diploid genotypes for grouped population samples
                    for base_name, haplotypes in sample_groups.items():
                        if 'h1' in haplotypes and 'h2' in haplotypes:
                            h1_record = haplotypes['h1']
                            h2_record = haplotypes['h2']
                            
                            # Extract genotypes and IDs
                            gt_h1_pop = h1_record.get('gt', '0')
                            gt_h2_pop = h2_record.get('gt', '0')
                            ids_h1_pop = h1_record.get('motif_ids_h', [])
                            ids_h2_pop = h2_record.get('motif_ids_h', [])
                            
                            # Compute diploid genotype
                            diploid_gt_pop = parsers.compute_diploid_genotype_assembly(gt_h1_pop, gt_h2_pop, ids_h1_pop, ids_h2_pop)
                            population_genotypes[base_name] = diploid_gt_pop
                        elif 'h1' in haplotypes:
                            population_genotypes[base_name] = haplotypes['h1'].get('gt', '0')
                        elif 'h2' in haplotypes:
                            population_genotypes[base_name] = haplotypes['h2'].get('gt', '0')
                    
                    if len(population_genotypes) > 1:
                        create_genotype_comparison_matrix(population_genotypes)
                        st.markdown("---")
                    
                    plots.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
                else:
                    st.info("no population data found")

            




def display_motif_legend(motifs, motif_colors, right_column):
        st.markdown("### Motif Legend")
        st.markdown('<div style="max-height:400px; overflow-y:scroll;">', unsafe_allow_html=True)  
        if  isinstance(motifs, tuple):
            motifs = list(motifs)
        elif not isinstance(motifs, list):
            motifs = [motifs]
        for idx, motif in enumerate(motifs):
            color = motif_colors[idx]
            st.markdown(
                f'<div id="legend-motif-{idx}" class="legend-item motif-{idx}" style="background-color:{color};color:white;padding:5px;margin-bottom:10px;border-radius:5px;'
                f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);'
                f' white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="{motif}">'
                f' Motif {idx}: {motif}</div>', unsafe_allow_html=True)
        st.markdown(
            f'<div id="legend-motif-interruption" class="legend-item motif-interruption" style="background-color:#FF0000;color:black;padding:5px;margin-bottom:10px;border-radius:5px;'
            f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
            f' Interruption</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)



def visulize_TR_with_dynamic_sequence(record, hgsvc_records, container, motif_colors, CN1_col, CN2_col, show_comparison):
        
        if record['motif_ids_h1'] == ['.'] and record['motif_ids_h2'] == ['.']:
            st.info("No motifs found in the region")
            return
        motif_names = record['motifs']
        reference_copy_number = record['ref_CN']
        motif_count_ref = parsers.count_motifs(record['motif_ids_ref'])
        found_motifs_ref = list(motif_count_ref.keys())
        found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]

        motif_count_h1 = parsers.count_motifs(record['motif_ids_h1'])
        found_motifs_h1 = list(motif_count_h1.keys())
        found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
        motif_count_h2 = parsers.count_motifs(record['motif_ids_h2'])
        found_motifs_h2 = list(motif_count_h2.keys())
        found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
        motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
        motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
        total_copy_number_h1 = str(record['spans'][1]).count('-')


    
 
        with container:
            st.markdown("""
                <style>
                    .stTabs [data-baseweb="tab-list"] {
                        gap: 6px;
                        background: #f8fafc;
                        padding: 4px 4px;
                        border-radius: 16px;
                        margin-bottom: 18px;
                        min-height: 22px;
                    }
                    .stTabs [data-baseweb="tab"] {
                        height: 28px;
                        min-height: 28px;
                        min-width: 38px;
                        white-space: pre-wrap;
                        background: white;
                        border-radius: 10px;
                        border: 2px solid #e2e8f0;
                        color: #64748b;
                        font-weight: 700;
                        font-size: 1.13rem;
                        padding: 0 13px;
                        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                        line-height: 1.33;
                    }
                    .stTabs [data-baseweb="tab"]:hover {
                        background: #f1f5f9;
                        border-color: #cbd5e1;
                        color: #475569;
                        transform: translateY(-1.8px) scale(1.04);
                        box-shadow: 0 2px 12px rgba(0,0,0,0.13);
                    }
                    .stTabs [aria-selected="true"] {
                        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                        color: white !important;
                        border-color: #667eea !important;
                        box-shadow: 0 5px 19px rgba(102, 126, 234, 0.23) !important;
                        transform: translateY(-2px) scale(1.045);
                        font-size: 1.23rem !important;
                    }
                    .stTabs [data-baseweb="tab-highlight"] {
                        background-color: transparent !important;
                    }
                </style>
            """, unsafe_allow_html=True)




        
            tab1, tab2, tab3 = st.tabs([
                "🧬 **Alleles**", 
                "🔄 **Alleles vs Ref**", 
                "🌍 **Alleles vs Pop**"
            ])
            with tab2:
                st.html('<div class="tab-content">')
                # Your tab1 content here
                motif_legend_html(record['motif_ids_ref'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("###  Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                display_dynamic_sequence_with_highlighted_motifs("Ref", record['ref_allele'], record['motif_ids_ref'], record['spans'][0], motif_colors, motif_names)
                alt_allele1 = record['alt_allele1']
                display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names, record['supported_reads_h1'])
                if record['alt_allele2'] != '':
                    alt_allele2 = record['alt_allele2'] 
                    display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names, record['supported_reads_h2'])
                st.html('</div>')
            with tab1:
                motif_legend_html(record['motif_ids_h1'], motif_colors, motif_names)
                
                # Add genotype visualization under motifs in region
                st.markdown("###  Genotype Information")
                display_genotype_card(record['gt'], "Current Sample", show_details=True)
                st.markdown("---")
                
                st.html('<div class="tab-content">')
                alt_allele1 = record['alt_allele1']
                display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names, record['supported_reads_h1'])

                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        with st.expander("Show motif bar plot for Allel1", expanded=False):
                            plot_motif_bar(motif_count_h1, motif_names, motif_colors)

                if record['alt_allele2'] != '':

                    alt_allele2 = record['alt_allele2'] 
                    display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names, record['supported_reads_h2'])
                    
                    plot_container_h2 = st.empty()
                    with plot_container_h2:
                        if show_comparison == False:
                            with st.expander("Show motif bar plot for Allel2", expanded=False):
                                plot_motif_bar(motif_count_h2, motif_names, motif_colors)
                st.html('</div>')
            with tab3:
                if hgsvc_records:
                    # Add genotype comparison for population data
                    # For current sample: compute diploid genotype from haplotypes (similar to cohort mode)
                    current_sample_gt = record['gt']
                    # Extract gt_h1 and gt_h2 from the GT field if available
                    # The GT field format is "0/1" or "0/0", so we need to split it
                    gt_parts = current_sample_gt.split('/') if '/' in str(current_sample_gt) else [current_sample_gt, current_sample_gt]
                    gt_h1 = gt_parts[0] if len(gt_parts) > 0 else "0"
                    gt_h2 = gt_parts[1] if len(gt_parts) > 1 else gt_parts[0]
                    
                    # Compute diploid genotype using the same logic as cohort mode
                    diploid_gt = parsers.compute_diploid_genotype_assembly(
                        gt_h1, gt_h2, 
                        record.get('motif_ids_h1', []), 
                        record.get('motif_ids_h2', [])
                    )
                    
                    population_genotypes = {"Current Sample": diploid_gt}
                    
                    # Handle population records - group by base name if they have _h1/_h2 suffixes
                    sample_groups = {}
                    for sample_name, sample_record in hgsvc_records.items():
                        if sample_name.endswith('_h1'):
                            base_name = sample_name[:-3]
                            if base_name not in sample_groups:
                                sample_groups[base_name] = {}
                            sample_groups[base_name]['h1'] = sample_record
                        elif sample_name.endswith('_h2'):
                            base_name = sample_name[:-3]
                            if base_name not in sample_groups:
                                sample_groups[base_name] = {}
                            sample_groups[base_name]['h2'] = sample_record
                        else:
                            # Sample doesn't have _h1/_h2 suffix, use as is
                            if 'gt' in sample_record:
                                population_genotypes[sample_name] = sample_record['gt']
                    
                    # Compute diploid genotypes for grouped population samples
                    for base_name, haplotypes in sample_groups.items():
                        if 'h1' in haplotypes and 'h2' in haplotypes:
                            h1_record = haplotypes['h1']
                            h2_record = haplotypes['h2']
                            
                            # Extract genotypes and IDs
                            gt_h1_pop = h1_record.get('gt', '0')
                            gt_h2_pop = h2_record.get('gt', '0')
                            ids_h1_pop = h1_record.get('motif_ids_h', [])
                            ids_h2_pop = h2_record.get('motif_ids_h', [])
                            
                            # Compute diploid genotype
                            diploid_gt_pop = parsers.compute_diploid_genotype_assembly(gt_h1_pop, gt_h2_pop, ids_h1_pop, ids_h2_pop)
                            population_genotypes[base_name] = diploid_gt_pop
                        elif 'h1' in haplotypes:
                            population_genotypes[base_name] = haplotypes['h1'].get('gt', '0')
                        elif 'h2' in haplotypes:
                            population_genotypes[base_name] = haplotypes['h2'].get('gt', '0')
                    
                    if len(population_genotypes) > 1:
                        create_genotype_comparison_matrix(population_genotypes)
                        st.markdown("---")
                    
                    plots.plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)
                else:
                    st.info("no population data found")

                st.html('</div>')
                

def visulize_cohort():
    if 'cohorts_records_map' in st.session_state:
        region_options = list(st.session_state.cohorts_records_map.values())
        
        # Safe index access with bounds checking
        regions_idx = st.session_state.get('regions_idx', 0)
        if regions_idx >= len(region_options):
            regions_idx = 0
            st.session_state.regions_idx = 0
        
        default_region = region_options[regions_idx] if region_options else ""
        
        # Cache the full options list to avoid recreating it
        if 'cached_region_options_cohort' not in st.session_state:
            st.session_state.cached_region_options_cohort = region_options
        
        # Single unified field: text input with autocomplete suggestions
        # Track if a selection was made
        if 'region_selected_cohort' not in st.session_state:
            st.session_state.region_selected_cohort = ""
        
        search_query = st.sidebar.text_input(
            "🔍 Search region:", 
            value=st.session_state.region_selected_cohort,
            key="region_search_cohort",
            help="Type to search and filter results",
            placeholder="Type to search..."
        )
        
        # Use markdown+CSS trick to make selectbox text black
        st.markdown("""
            <style>
            /* Ensure all selectbox texts are black, regardless of state */
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] span {
                color: black !important;
            }
            [data-testid="stSidebar"] .stSelectbox label, 
            [data-testid="stSidebar"] .stSelectbox div[role="listbox"] span {
                color: black !important;
            }
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] input,
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="combobox"] span,
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="button"] span,
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="option"] span,
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[role="listbox"] span {
                color: black !important;
            }
            /* Also set all the selectbox selected value text to black */
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] input {
                color: black !important;
                font-weight: 700;
            }
            /* For v1.31+ (possible dark mode or material theme changes): */
            [data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] div[aria-selected="true"] span {
                color: black !important;
            }
            </style>
        """, unsafe_allow_html=True)

        
        # Filter and show suggestions only when typing
        if search_query:
            search_lower = search_query.lower()
            filtered = [r for r in st.session_state.cached_region_options_cohort if search_lower in r.lower()]
            filtered_regions = filtered[:10]  # Limit to 10 results
            total_matches = len(filtered)
            
            # Show filtered suggestions as a dropdown without label to feel unified
            if filtered_regions:
                region = st.sidebar.selectbox(
                    " ",  # Empty label
                    filtered_regions, 
                    index=0,
                    key="region_suggest",
                    help="Select from suggestions",
                    label_visibility="collapsed"
                )
                # Update the text input with the selected region
                st.session_state.region_selected_cohort = region
            else:
                region = search_query
                
            # Show match count
            if total_matches > 10:
                st.sidebar.markdown(f"<span style='font-size:11px; color:orange;'>Showing 10 of {total_matches:,} matches</span>", unsafe_allow_html=True)
            else:
                st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>{total_matches:,} matches</span>", unsafe_allow_html=True)
        else:
            # When not typing, use the current region from session state
            region = default_region
            total_matches = len(st.session_state.cached_region_options_cohort)
            st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>Search to find from {total_matches:,} regions</span>", unsafe_allow_html=True)

        if 'regions_idx' not in st.session_state:
            st.session_state.regions_idx = 0
        
        # Beautiful navigation buttons in sidebar
        st.sidebar.markdown("""
            <div style="
                display: flex;
                align-items: center;
                gap: 10px;
                background: linear-gradient(92deg, #667eea 0%, #a084e8 100%);
                padding: 10px 18px;
                border-radius: 14px;
                box-shadow: 0 2px 14px rgba(102, 126, 234, 0.10);
                margin-bottom: 8px;
                ">
                <span style="
                    font-size: 28px;
                    color: #4338ca;
                    margin-right: 6px;
                    filter: drop-shadow(0 2px 6px rgba(65,0,140,0.18));
                    ">🧭</span>
                <span style="
                    font-size: 1.25rem; 
                    font-weight: 700; 
                    letter-spacing: 0.03em; 
                    color: #312e81;
                    ">Region Navigation</span>
            </div>
        """, unsafe_allow_html=True)
        nav_col1, nav_col2 = st.sidebar.columns(2, gap="small")
        with nav_col1:
            if st.button("◀ Previous", use_container_width=True, key="prev_region"):
                region = None
                st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)
        with nav_col2:
            if st.button("Next ▶", use_container_width=True, key="next_region"):
                region = None
                st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len( st.session_state.cohorts_records_map )-1)
        
        # Add beautiful styling for sidebar navigation buttons
        st.markdown("""
            <style>
                /* Beautiful sidebar navigation buttons */
                [data-testid="stSidebar"] button[kind="secondary"] {
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                    color: white !important;
                    border: none !important;
                    border-radius: 12px !important;
                    font-weight: 700 !important;
                    font-size: 15px !important;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
                    box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3) !important;
                }
                [data-testid="stSidebar"] button[kind="secondary"]:hover {
                    transform: translateY(-2px) scale(1.05) !important;
                    box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4) !important;
                    background: linear-gradient(135deg, #764ba2 0%, #667eea 100%) !important;
                }
            </style>
            <script>
            (function() {
                function styleNavButtons() {
                    const sidebarButtons = document.querySelectorAll('[data-testid="stSidebar"] button');
                    sidebarButtons.forEach(btn => {
                        const text = btn.textContent || btn.innerText || '';
                        if (text.includes('◀ Previous') || text.includes('Next ▶')) {
                            btn.style.background = 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)';
                            btn.style.color = 'white';
                            btn.style.border = 'none';
                            btn.style.borderRadius = '12px';
                            btn.style.fontWeight = '700';
                            btn.style.fontSize = '15px';
                            btn.style.transition = 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)';
                            btn.style.boxShadow = '0 4px 15px rgba(102, 126, 234, 0.3)';
                        }
                    });
                }
                styleNavButtons();
                setTimeout(styleNavButtons, 100);
                setTimeout(styleNavButtons, 500);
                const observer = new MutationObserver(styleNavButtons);
                observer.observe(document.body, { childList: true, subtree: true });
            })();
            </script>
        """, unsafe_allow_html=True)
        if region and region != st.session_state.get('previous_region', None):
            try:
                chr_input, start_end_input = region.split(':')
                start_input, end_input = map(int, start_end_input.split('-'))

                region = f"{chr_input}:{start_input}-{end_input}"
                st.session_state.regions_idx = list(st.session_state.cohorts_records_map.values()).index(region)
            except:
                try:
                    chr_input, start_input, end_input = re.split(r'\s+', region)
                    start_input, end_input = int(start_input), int(end_input)
                    region = f"{chr_input}:{start_input}-{end_input}"
                except:
                    st.sidebar.info("Invalid region format, showing the first record")
                    region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
                
            st.session_state.previous_region = region
        else:
            #st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len( st.session_state.cohorts_records_map )-1)
            region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
        mode_placeholder = st.empty()
        
        # Create container for region display
        middel, _ = st.columns([1, 0.1], gap="small")
        render_region_display(middel, region)
        st.markdown("""
                <style>
                :root {
                    --region-color-light: black;
                    --region-color-dark: white;
                }
                /* Default style for light mode */
                .region-container {
                    color: var(--region-color-light);
                }
                /* Apply different color for dark mode */
                @media (prefers-color-scheme: dark) {
                    .region-container {
                        color: var(--region-color);
                    }
                }
                </style>
            """, unsafe_allow_html=True)

        st.session_state.cohort_results = utils.get_results_cohort(
            region, st.session_state.cohort_files, st.session_state.cohort_file_paths, 
            st.session_state.cohort_mode,
            parsers.parse_record,
            parsers.parse_record_assembly
        )
        if 'cohort_results' in st.session_state:   
            region = st.session_state.regions_idx
            if (st.session_state.cohort_results == {}):
                st.warning(f"No motifs found in the region: {region}")
                st.stop()
            else:
                    plots.plot_Cohort_results(st.session_state.cohort_results, st.session_state.cohort_mode)
        else:
            st.stop()



def visulize_region():
    if 'regions_idx' not in st.session_state:
        st.session_state.regions_idx = 0
    region_options = list(st.session_state.records_map.values())
    
    # Safe index access with bounds checking
    regions_idx = st.session_state.get('regions_idx', 0)
    if regions_idx >= len(region_options):
        regions_idx = 0
        st.session_state.regions_idx = 0
    
    default_region = region_options[regions_idx] if region_options else ""
    st.sidebar.markdown("### Select Region to Visualize")
    
    # Cache the full options list to avoid recreating it
    if 'cached_region_options' not in st.session_state:
        st.session_state.cached_region_options = region_options
    
    # Track if a selection was made
    if 'region_selected_ind' not in st.session_state:
        st.session_state.region_selected_ind = ""
    
    # Searchable field: text input that filters selectbox
    search_query = st.sidebar.text_input(
        "🔍 Search region:", 
        value=st.session_state.region_selected_ind,
        key="region_search",
        help="Type to search and filter results",
        placeholder="Type to search..."
    )
    
    # Filter options based on search and show only 10 results
    if search_query:
        search_lower = search_query.lower()
        filtered = [r for r in st.session_state.cached_region_options if search_lower in r.lower()]
        filtered_regions = filtered[:10]  # Show only 10 results
        total_matches = len(filtered)
        
        # Show the filtered selectbox without label to feel unified
        if filtered_regions:
            region = st.sidebar.selectbox(
                " ",  # Empty label
                filtered_regions, 
                index=0,
                key="region_suggest",
                help="Select from filtered results",
                label_visibility="collapsed"
            )
            # Update the text input with the selected region
            st.session_state.region_selected_ind = region
        else:
            region = search_query
            
        if total_matches > 10:
            st.sidebar.markdown(f"<span style='font-size:11px; color:orange;'>Showing 10 of {total_matches:,} matches</span>", unsafe_allow_html=True)
        else:
            st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>{total_matches:,} matches</span>", unsafe_allow_html=True)
    else:
        # No search - show default region
        region = default_region
        total_matches = len(st.session_state.cached_region_options)
        st.sidebar.markdown(f"<span style='font-size:11px; color:white;'>Search to find from {total_matches:,} regions</span>", unsafe_allow_html=True)
    
    display_option = st.sidebar.radio("Select Display Type", 
                                        ("Sequence with Highlighted Motifs", "Bars"))

    # Beautiful navigation buttons in sidebar
    st.sidebar.markdown("### 🧭 Navigation")
    nav_col1, nav_col2 = st.sidebar.columns(2, gap="small")
    with nav_col1:
        if st.button("◀ Previous", use_container_width=True, key="prev_individual"):
            region = None
            st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

    with nav_col2:
        if st.button("Next ▶", use_container_width=True, key="next_individual"):
            region = None
            st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.records_map) - 1)
    
    # Add beautiful styling for sidebar navigation buttons
    st.markdown("""
        <style>
            /* Beautiful sidebar navigation buttons */
            [data-testid="stSidebar"] button[kind="secondary"] {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
                color: white !important;
                border: none !important;
                border-radius: 12px !important;
                font-weight: 700 !important;
                font-size: 15px !important;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
                box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3) !important;
            }
            [data-testid="stSidebar"] button[kind="secondary"]:hover {
                transform: translateY(-2px) scale(1.05) !important;
                box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4) !important;
                background: linear-gradient(135deg, #764ba2 0%, #667eea 100%) !important;
            }
        </style>
        <script>
        (function() {
            function styleNavButtons() {
                const sidebarButtons = document.querySelectorAll('[data-testid="stSidebar"] button');
                sidebarButtons.forEach(btn => {
                    const text = btn.textContent || btn.innerText || '';
                    if (text.includes('◀ Previous') || text.includes('Next ▶')) {
                        btn.style.background = 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)';
                        btn.style.color = 'white';
                        btn.style.border = 'none';
                        btn.style.borderRadius = '12px';
                        btn.style.fontWeight = '700';
                        btn.style.fontSize = '15px';
                        btn.style.transition = 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)';
                        btn.style.boxShadow = '0 4px 15px rgba(102, 126, 234, 0.3)';
                    }
                });
            }
            styleNavButtons();
            setTimeout(styleNavButtons, 100);
            setTimeout(styleNavButtons, 500);
            const observer = new MutationObserver(styleNavButtons);
            observer.observe(document.body, { childList: true, subtree: true });
        })();
        </script>
    """, unsafe_allow_html=True)
    
    middel, spacer = st.columns([1, 0.1], gap="small")
    REF, CN1_col, CN2_col = st.columns([1, 1, 1])
    if region and region != st.session_state.get('previous_region', None):
        try:
            chr_input, start_end_input = region.split(':')
            start_input, end_input = map(int, start_end_input.split('-'))
        

            input_region = f"{chr_input}:{start_input}-{end_input}"
            record_key = st.session_state.records[input_region]
            st.session_state.regions_idx = list(st.session_state.records_map.values()).index(input_region)
                

        except:
            try:
                chr_input, start_input, end_input = re.split(r'\s+', region)
                start_input, end_input = int(start_input), int(end_input)
                input_region = f"{chr_input}:{start_input}-{end_input}"
                record_key = st.session_state.records[input_region]
                st.session_state.regions_idx = list(st.session_state.records_map.values()).index(input_region)
            except:
                st.sidebar.info("Invalid region format, showing the first record")
                record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]
    else:
        record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]
    st.session_state.previous_region = region
    record = parsers.parse_record(st.session_state.vcf_file_path, record_key)
    if record is None:
        st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
        st.stop()
    hgsvc_records = utils.get_results_hgsvc_pop(record_key, st.session_state.files, st.session_state.file_paths, parsers.parse_record_assembly)
    if len(record["motif_ids_h1"]) == 0 and len(record["motif_ids_h2"]) == 0:
        st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
        st.stop()

    render_region_display(middel, st.session_state.records_map[st.session_state.regions_idx])


    st.markdown("""
                    <style>
                    :root {
                        --region-color-light: black;
                        --region-color-dark: white;
                    }
                    /* Default style for light mode */
                    .region-container {
                        color: var(--region-color-light);
                    }
                    /* Apply different color for dark mode */
                    @media (prefers-color-scheme: dark) {
                        .region-container {
                            color: var(--region-color);
                        }
                    }
                    </style>
                """, unsafe_allow_html=True)

    container = st.container()
    motif_colors = utils.get_color_palette(len(record['motifs']))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
        


    if display_option == "Sequence with Highlighted Motifs":
        visulize_TR_with_dynamic_sequence(record, hgsvc_records, container, motif_colors, CN1_col, CN2_col, st.session_state.get('show_comparison', False))

    elif display_option == "Bars":
        display_motifs_with_bars(record, container, motif_colors, CN1_col, CN2_col, st.session_state.get('show_comparison', False), hgsvc_records)



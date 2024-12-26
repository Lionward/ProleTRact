import streamlit as st
import pysam
import re
import matplotlib.pyplot as plt
import pandas as pd
import altair as alt
import os
import pysam
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(layout="wide")


st.session_state.pathogenic_TRs = {
    "chr1:146228800-146228821": {"gene": "NOTCH2NLA", "pathogenicity_threshold": [61]},
    "chr1:149390802-149390841": {"gene": "NOTCH2NL", "pathogenicity_threshold": [66, 517]},
    "chr10:79826383-79826404": {"gene": "NUTM2B-AS1", "pathogenicity_threshold": [16, 160]},
    "chr10:93702522-93702547": {"gene": "FRA10AC1", "pathogenicity_threshold": [0]},
    "chr11:66744821-66744850": {"gene": "C11ORF80", "pathogenicity_threshold": [0]},
    "chr11:119206289-119206322": {"gene": "CBL", "pathogenicity_threshold": [101]},
    "chr12:6936716-6936773": {"gene": "ATN1", "pathogenicity_threshold": [49, 93]},
    "chr12:50505001-50505022": {"gene": "DIP2B", "pathogenicity_threshold": [151]},
    "chr12:111598949-111599018": {"gene": "ATXN2", "pathogenicity_threshold": [33, 200]},
    "chr13:70139353-70139428": {"gene": "ATXN8", "pathogenicity_threshold": [74, 1300]},
    "chr13:99985448-99985494": {"gene": "ZIC2", "pathogenicity_threshold": [25]},
    "chr14:23321472-23321490": {"gene": "PABPN1", "pathogenicity_threshold": [18]},
    "chr14:92071009-92071042": {"gene": "ATXN3", "pathogenicity_threshold": [53, 87]},
    "chr15:22786677-22786701": {"gene": "NIPA1", "pathogenicity_threshold": [9]},
    "chr16:17470907-17470923": {"gene": "XYLT1", "pathogenicity_threshold": [120]},
    "chr16:24613439-24613530": {"gene": "TNRC6A", "pathogenicity_threshold": [28]},
    "chr16:66490398-66490467": {"gene": "BEAN1", "pathogenicity_threshold": [500, 760]},
    "chr16:87604287-87604329": {"gene": "JPH3", "pathogenicity_threshold": [40, 59]},
    "chr18:55586155-55586227": {"gene": "TCF4", "pathogenicity_threshold": [51]},
    "chr19:13207858-13207897": {"gene": "CACNA1A", "pathogenicity_threshold": [19, 33]},
    "chr19:14496041-14496074": {"gene": "GIPC1", "pathogenicity_threshold": [70, 164]},
    "chr19:18786034-18786050": {"gene": "COMP", "pathogenicity_threshold": [6]},
    "chr19:45770204-45770264": {"gene": "DMPK", "pathogenicity_threshold": [50, 10000]},
    "chr2:96197066-96197122": {"gene": "STARD7", "pathogenicity_threshold": [150, 460]},
    "chr2:100104799-100104824": {"gene": "AFF3", "pathogenicity_threshold": [0]},
    "chr2:176093058-176093104": {"gene": "HOXD13", "pathogenicity_threshold": [22]},
    "chr2:190880872-190880920": {"gene": "GLS", "pathogenicity_threshold": [680]},
    "chr20:2652733-2652775": {"gene": "NOP56", "pathogenicity_threshold": [650, 2500]},
    "chr21:43776443-43776479": {"gene": "CSTB", "pathogenicity_threshold": [30, 125]},
    "chr22:45795354-45795424": {"gene": "ATXN10", "pathogenicity_threshold": [280, 4500]},
    "chr3:63912684-63912726": {"gene": "ATXN7", "pathogenicity_threshold": [34, 460]},
    "chr3:129172576-129172732": {"gene": "CNBP", "pathogenicity_threshold": [50, 11000]},
    "chr3:138946020-138946063": {"gene": "FOXL2", "pathogenicity_threshold": [22]},
    "chr3:183712187-183712223": {"gene": "YEATS2", "pathogenicity_threshold": [201]},
    "chr4:3074876-3074966": {"gene": "HTT", "pathogenicity_threshold": [36, 250]},
    "chr4:39348424-39348479": {"gene": "RFC1", "pathogenicity_threshold": [400, 2000]},
    "chr4:41745972-41746032": {"gene": "PHOX2B", "pathogenicity_threshold": [24]},
    "chr4:159342526-159342617": {"gene": "RAPGEF2", "pathogenicity_threshold": [28]},
    "chr5:10356346-10356412": {"gene": "MARCHF6", "pathogenicity_threshold": [700, 1035]},
    "chr5:146878727-146878757": {"gene": "PPP2R2B", "pathogenicity_threshold": [51, 78]},
    "chr6:16327633-16327723": {"gene": "ATXN1", "pathogenicity_threshold": [39, 91]},
    "chr6:45422750-45422802": {"gene": "RUNX2", "pathogenicity_threshold": [27]},
    "chr6:170561906-170562017": {"gene": "TBP", "pathogenicity_threshold": [43, 66]},
    "chr7:27199825-27199862": {"gene": "HOXA13", "pathogenicity_threshold": [24]},
    "chr7:55887601-55887640": {"gene": "ZNF713", "pathogenicity_threshold": [None]},
    "chr8:104588972-104589000": {"gene": "LRP12", "pathogenicity_threshold": [90, 130]},
    "chr8:118366815-118366919": {"gene": "SAMD12", "pathogenicity_threshold": [105, 3680]},
    "chr9:27573528-27573546": {"gene": "C9ORF72", "pathogenicity_threshold": [24, 4000]},
    "chr9:69037261-69037304": {"gene": "FXN", "pathogenicity_threshold": [66, 1300]},
    "chrX:25013649-25013698": {"gene": "ARX", "pathogenicity_threshold": [17, 27]},
    "chrX:67545316-67545385": {"gene": "AR", "pathogenicity_threshold": [38, 68]},
    "chrX:140504316-140504362": {"gene": "SOX3", "pathogenicity_threshold": [15]},
    "chrX:147912050-147912110": {"gene": "FMR1", "pathogenicity_threshold": [55, 200, 3000]},
    "chrX:148500631-148500691": {"gene": "AFF2", "pathogenicity_threshold": [200]},
    "chrX:149631763-149631782": {"gene": "TMEM185A", "pathogenicity_threshold": [0]},
}



# JavaScript code to check the color scheme
# ms = st.session_state
# if "themes" not in ms: 
#   ms.themes = {"current_theme": "light",
#                     "refreshed": True,
                    
#                     "light": {"theme.base": "dark",
#                               "theme.backgroundColor": "black",
#                               "theme.primaryColor": "#c98bdb",
#                               "theme.secondaryBackgroundColor": "#5591f5",
#                               "theme.textColor": "white",
#                               "theme.textColor": "white",
#                               "button_face": "🌜"},

#                     "dark":  {"theme.base": "light",
#                               "theme.backgroundColor": "white",
#                               "theme.primaryColor": "#5591f5",
#                               "theme.secondaryBackgroundColor": "#82E1D7",
#                               "theme.textColor": "#0a1464",
#                               "button_face": "🌞"},
#                     }
  

# def ChangeTheme():
#   st.write("Changing theme")
#   st.write("previous theme", ms.themes["current_theme"])

#   previous_theme = ms.themes["current_theme"]
#   tdict = ms.themes["light"] if ms.themes["current_theme"] == "light" else ms.themes["dark"]
#   for vkey, vval in tdict.items(): 
#     if vkey.startswith("theme"): st._config.set_option(vkey, vval)

#   ms.themes["refreshed"] = False
#   if previous_theme == "dark": 
#     ms.themes["current_theme"] = "light"
#   elif previous_theme == "light":
#     ms.themes["current_theme"] = "dark"
#     st.write("current theme", ms.themes["current_theme"])
    
# btn_face = ms.themes["light"]["button_face"] if ms.themes["current_theme"] == "light" else ms.themes["dark"]["button_face"]
# st.button(btn_face, on_click=ChangeTheme)

# if ms.themes["refreshed"] == False:
#     ms.themes["refreshed"] = True
#     st.rerun()

@st.cache_data()
def parse_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    records_ids = {}
    records_map = {}
    idx = 0
    for  rec in vcf.fetch():
        records_ids[rec.id] = f"{rec.chrom}:{rec.pos}-{rec.stop}"
        records_map[idx] = rec.id
        idx += 1
        
    return records_ids, records_map

def load_vcf(vcf_file):
    # check if the vcf file has an index file
    
    return pysam.VariantFile(vcf_file)

def parse_record_assembly(vcf,region):        
        chr,start_end = region.split(":")
        start,end = start_end.split("-")
        start = int(start) -1
        end = int(end) -1
        region = f"{chr}:{start}-{end}"
        
        for rec in vcf.fetch(region=region):
            break
        try:
            ids_h = rec.info['MOTIF_IDs_H']
            ref_CN = rec.info['CN_ref']
            alt_allele = "."

            ref_allele = rec.ref
            if rec.alts != None:
                if len(rec.alts) > 0:
                    alt_allele = rec.alts[0]
                    if alt_allele == '.':
                        alt_allele = ''
                else:
                    alt_allele = '' 
            if ids_h is None:
                ids_h = []
        
            CN_H = rec.info['CN_hap']
    
            motif_names = rec.info['MOTIFS']
            if isinstance(motif_names, tuple):
                motif_names = list(motif_names)
            elif not isinstance(motif_names, list):
                motif_names = [motif_names]
        
        
            record = {
                    'chr': rec.chrom,
                    'pos': rec.pos,
                    'stop': rec.stop,
                    'motifs': motif_names,
                    'motif_ids_h': ids_h,
                    'motif_ids_ref': rec.info['MOTIF_IDs_REF'],
                    'ref_CN': ref_CN,
                    'CN_H': CN_H,
                    'spans': rec.samples[0]['SP'],
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele,
                }
        except:
            record = {
                    'chr': "",
                    'pos': -1,
                    'stop': -1,
                    'motifs': [],
                    'motif_ids_h': [],
                    'motif_ids_ref': [],
                    'ref_CN': 0,
                    'CN_H': 0,
                    'spans': [],
                    'ref_allele': '',
                    'alt_allele': '',
                }
        return record


def get_results_hgsvc_pop(region, files, file_paths):
    samples_results = {}
    for i in range(len(files)):
        sample_name = file_paths[i].split(".")[0]

        record = parse_record_assembly(st.session_state.files[i], region)
        samples_results[sample_name] = record
    return samples_results

def get_results_cohort(region, files, file_paths):
    samples_results = {}
    chrom, start_end = region.split(":")
    start, end = start_end.split("-")
    start = int(start) - 1
    end = int(end) - 1
    region = f"{chrom}:{start}-{end}"

    for i in range(len(files)):
        sample_name = file_paths[i].split(".")[0]
        record = parse_record(files[i], region)
        samples_results[sample_name] = record
    return samples_results

def parse_record(vcf_file,region):
    is_vcf_file_str = isinstance(vcf_file, str)
    if is_vcf_file_str:
        vcf = pysam.VariantFile(vcf_file)
    else:
        vcf = vcf_file
    rec = vcf.fetch(region=region)
    # get the record with the id
    for rec in vcf.fetch(region=region):
        break
    if st.session_state.show_vcf:
        st.write(rec)
    ids_h1 = rec.info['MOTIF_IDs_H1']
    ids_h2 = rec.info['MOTIF_IDs_H2']
    ref_CN = rec.info['CN_ref']
    alt_allele1 = "."
    alt_allele2 = "."
    ref_allele = rec.ref
    if rec.alts != None:
        if len(rec.alts) > 0:
            alt_allele1 = rec.alts[0]
            if alt_allele1 == '.':
                alt_allele1 = ''
        if len(rec.alts) > 1:
            alt_allele2 = rec.alts[1]
        else:
            if ids_h1 == ids_h2:
                alt_allele2 = alt_allele1
            if rec.samples[0]["DP"][1] == 0:
                alt_allele2 = alt_allele1
                ids_h2 = ids_h1
            #alt_allele2 = ''
    if not ids_h1:
        ids_h1 = []
    if not ids_h2:
        ids_h2 = []
   
    cn_h1 = rec.info['CN_H1']
    cn_h2 = rec.info['CN_H2']

    if cn_h1 == cn_h2 and len(ids_h2) == 0:
        spans = [ rec.samples[0]['SP'][0], rec.samples[0]['SP'][1], rec.samples[0]['SP'][1]]
        ids_h2 = ids_h1
    else:
        spans = [ rec.samples[0]['SP'][0], rec.samples[0]['SP'][1], rec.samples[0]['SP'][2]]
    # if alt_allele2 == 'N':
    #     spans = [ rec.samples[0]['SP'][0], rec.samples[0]['SP'][1]]
    #     alt_allele2 = ""

    motif_names = rec.info['MOTIFS']
    if isinstance(motif_names, tuple):
        motif_names = list(motif_names)
    elif not isinstance(motif_names, list):
        motif_names = [motif_names]

    if len(rec.samples[0]['SP']) >2 :
        spans= [ rec.samples[0]['SP'][0], rec.samples[0]['SP'][1], rec.samples[0]['SP'][2]]
    elif (len(rec.samples[0]['SP']) == 2 and alt_allele2 == alt_allele1):
        spans = [ rec.samples[0]['SP'][0], rec.samples[0]['SP'][1], rec.samples[0]['SP'][1]]
    elif (len(rec.samples[0 ]['SP']) == 2 and alt_allele2 != alt_allele1):
        spans = [ rec.samples[0]['SP'][0], rec.samples[0]['SP'][1], ""]
    if ids_h1 == ids_h2:
        spans[2] = spans[1]
        alt_allele2 = alt_allele1
    if alt_allele2 == '':
        spans[2] = ""
    if alt_allele1 == '':
        spans[1] = ""

    record = {
            'chr': rec.chrom,
            'pos': rec.pos,
            'stop': rec.stop,
            'motifs': motif_names,
            'motif_ids_h1': ids_h1,
            'motif_ids_h2': ids_h2,
            'motif_ids_ref': rec.info['MOTIF_IDs_REF'],
            'ref_CN': ref_CN,
            'CN_H1': cn_h1,
            'CN_H2': cn_h2,
            'spans': spans,
            'ref_allele': ref_allele,
            'alt_allele1': alt_allele1,
            'alt_allele2': alt_allele2
        }
    
    return record

def parse_motif_range(motif_range):
    if motif_range is {}:
        return []
    pattern = re.compile(r'\((\d+)-(\d+)\)')
    matches = pattern.findall(motif_range)
    ranges = [(int(start)-1, int(end)-1) for start, end in matches]
    return ranges

# Function to generate a list of n visually distinct colors using a color map
def get_color_palette(n):
    
    #cmap = plt.get_cmap('tab20') # more options 'viridis', 'plasma', 'Set1', etc.
    # add 100 more colors if n > 20
    colors_names = ['tab20', 'tab20c', 'tab20b', 'tab20', 'tab20c', 'tab20', 'tab20c', 'gist_rainbow']
    a = 0
    cmap_colors = []
    if n > 20:
        while (len(cmap_colors) < n):
            camp = plt.get_cmap(colors_names[a])
            cmap_colors += [camp(i) for i in range(camp.N)]
            a += 1
    else:
        camp = plt.get_cmap('tab20')
        cmap_colors = [camp(i) for i in range(n)]
    
    cmap_colors = cmap_colors[:n]
    # colors = [cmap(i) for i in range(n)]
    # st.write(len(colors))

    return ['#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in cmap_colors]

def display_dynamic_sequence_with_highlighted_motifs(sequence_name, sequence, motif_ids, spans, motif_colors, motif_names):
    ranges = parse_motif_range(spans)
    highlighted_sequence = ""
    previous_end = 0
   
    # Loop through each motif in the sequence and highlight
    for idx, (start, end) in enumerate(ranges):
        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        motif_name = motif_names[int(motif)]
        # Add the sequence before the motif and highlight interruptions in red
        if start > previous_end:
            interruption_sequence = sequence[previous_end:start]
            highlighted_sequence += (
                f"<span style='background-color:#FF0000; padding:2px; border-radius:4px;' title='Interruption'>"
                f"{interruption_sequence}</span>"
            )
        
        # Highlight the motif in the sequence
        motif_sequence = sequence[start:end+1]
        highlighted_sequence += (
            f"<span style='background-color:{color}; padding:2px; border-radius:4px;' title='Motif: {motif_name}'>"
            f"{motif_sequence}</span>"
        )

        previous_end = end + 1

    # Add remaining sequence after the last motif and highlight interruptions in red
    if previous_end < len(sequence):
        interruption_sequence = sequence[previous_end:]
        highlighted_sequence += (
            f"<span style='background-color:#FF0000; padding:2px; border-radius:4px;'title='Interruption'>"
            f"{interruption_sequence}</span>"
        )
    if sequence_name == "Ref":
        sequence_name += "seq"
    # Add scrollable container for the sequence display
    st.markdown(f"""
        <div style="display: flex; align-items: right;">
            <div style="font-family:monospace; font-size:16px; padding:10px; border:1px solid black; border-radius:8px; margin-right: 10px; text-align: right;">
                <strong>{sequence_name}</strong>
            </div>
            <div style="font-family:monospace; font-size:16px; width:100%; max-height:120px; overflow-x:auto; white-space:nowrap; padding:10px; border:1px solid black; border-radius:8px;">
                {highlighted_sequence}
            </div>
        </div>
    """, unsafe_allow_html=True)
# Function to display motifs relative to the sequence, with spans aligned horizontally
def display_motifs_with_bars(record, left_column, right_column,motif_colors,CN1_col,CN2_col, show_comparison):
    motif_names, motif_count_h1, motif_count_h2 = parse_motif_in_region(record)
    # iterate over the motifs and set them to 0 if they are not in the motif_count
    CN1_col.markdown(f""" 
        <div style="font-size: 20px; color: #FF5733;">
            <strong>Allele 1 Total copy number:</strong> {str(record['spans'][1]).count('-')}
        </div>
    """, unsafe_allow_html=True)

    with right_column:
        
        display_motif_legend(motif_names, motif_colors, right_column)

    if record['alt_allele2'] != '' and record['spans'][2] != None:
        CN2_col.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 2 Total copy number:</strong> {str(record['spans'][2]).count('-')}
            </div>
        """, unsafe_allow_html=True)
    else:
        CN2_col.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 2 Total copy number:</strong> &lt;DEL&gt;
            </div>
        """, unsafe_allow_html=True)


    with left_column:
        tab1, tab2 , tab3 = st.tabs(["Alleles", "Alleles vs Ref", "Alleles vs Pop"])
        # make the tabs font size bigger
        st.markdown(
            "<style>.tab-content {font-size: 20px;}</style>",
            unsafe_allow_html=True,
        )
        with tab2:
            display_motifs_as_bars("Ref", motif_colors, record['motif_ids_ref'], record['spans'][0], record['ref_allele'], motif_names)
            if record['alt_allele1'] != '' and record['spans'][1] != None:
                display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names)
            if record['alt_allele2'] != '' and record['spans'][2] != None:
                display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names)
        with tab1:
            if record['alt_allele1'] != '' and record['spans'][1] != None:
                display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names)
                plot_container_h1 = st.empty()
                with plot_container_h1:
                    if show_comparison == False:
                        plot_motif_bar(motif_count_h1, motif_names, motif_colors)
            
            if record['alt_allele2'] != '' and record['spans'][2] != None:
                display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names)
                plot_container_h2 = st.empty()

                with plot_container_h2:
                    if show_comparison == False:
                        plot_motif_bar(motif_count_h2, motif_names, motif_colors)

        with tab3:
            # add the alleles 
            plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)

def parse_motif_in_region(record):
    motif_names = record['motifs']
    motif_count_h1 = {}
    motif_count_h2 = {}

    if record['alt_allele2'] != '.' and  record['spans'][1] !=None:
        motif_count_ref = count_motifs(record['motif_ids_ref'], record['spans'][0])
        found_motifs_ref = list(motif_count_ref.keys())
        found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]
        motif_count_h1 = count_motifs(record['motif_ids_h1'], record['spans'][1])
        found_motifs_h1 = list(motif_count_h1.keys())
        found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
        motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
        motif_count_h2 = {}
    if record['alt_allele2'] != '.' and  record['spans'][2] != None:
        motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][2])
        found_motifs_h2 = list(motif_count_h2.keys())
        found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
        motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
    return motif_names,motif_count_h1,motif_count_h2




def display_motifs_as_bars(sequence_name, motif_colors, motif_ids, spans, sequence, motif_names):
    sequence_length = len(sequence)
    ranges = parse_motif_range(spans)

    # Ensure motif_names is a list
    if not isinstance(motif_names, list):
        motif_names = [motif_names]
    
    # Initialize bar container for visualizing motifs
    bar_container = "<div style='width:100%; position: relative; height: 30px; border:2px solid black; border-radius: 8px;'>"
    previous_end = 0
    gap = 0.05  # Small gap between motifs for better visibility

    # Loop through each motif in the sequence
    for idx, (start, end) in enumerate(ranges):
        # st.sidebar.write(len(sequence))
        # st.sidebar.write(len(motif_ids))

        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        span_length = end - start + 1  # Calculate length of motif span

        # Ensure the motif fits within the sequence length
        if start >= 0 and end <= sequence_length:
            # Handle interruptions between motifs
            if start > previous_end:
                interruption_width = (start - previous_end) / sequence_length * 100
                interruption_start = previous_end / sequence_length * 100
                # Make it also hoverable
                bar_container += (
                    f"<div style='position:absolute; background-color:#FF0000; left:{interruption_start}%; "
                    f"width:{interruption_width}%; height:28px; top:-1px; border-radius:6px; border:1px solid black; "
                    f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2); cursor: pointer;' "
                    f"title='Interruption: {sequence[previous_end:start]}'>"
                    f"</div>"
                )

            # Calculate relative width and position for the motif bar
            relative_width = (span_length / sequence_length) * 100 - gap
            relative_start = (start / sequence_length) * 100

            # Add each motif bar to the container with hover and color effects
            bar_container += (
                f"<div class='hoverable-div' style='position:absolute; background-color:{color}; left:{relative_start}%; "
                f"width:{relative_width}%; height:28px; top:-1px; border-radius:6px; border:1px solid black; "
                f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2); cursor: pointer;' "
                f"title='Motif: {motif_names[int(motif)]}'>"
                f"</div>"
            )

            previous_end = end + 1

    # Add interruption bar after the last motif, if any
    if previous_end < sequence_length:
        interruption_width = (sequence_length - previous_end) / sequence_length * 100
        interruption_start = previous_end / sequence_length * 100
        bar_container += (
            f"<div style='position:absolute; background-color:#FF0000; left:{interruption_start}%; "
            f"width:{interruption_width}%; height:28px; top:-1px; border-radius:6px; border:1px solid black; "
            f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2); cursor: pointer;' "
            f"title='Interruption: {sequence[previous_end:]}'>"
            f"</div>"
        )

    # Close the div container
    bar_container += "</div>"
    if sequence_name == "Ref":
        sequence_name += "seq"
    # Render the motif bars in Streamlit with the sequence name on the left
    st.markdown(f"""
        <div style="display: flex; align-items: center;">
            <div style="font-family:monospace; font-size:16px; padding:10px; border:1px solid black; border-radius:8px; margin-right: 10px;">
                <strong>{sequence_name}</strong>
            </div>
            {bar_container}
        </div>
    """, unsafe_allow_html=True)

def plot_motif_bar(motif_count, motif_names, motif_colors=None):
    # Convert motif indices to names and store counts
    motif_labels = []
    motif_counts = []
    for label, value in sorted(motif_count.items()):
        motif_name = motif_names[int(label)]
        if motif_name:
            motif_labels.append(motif_name)
            motif_counts.append(value)
        
    # Create DataFrame with motif labels and counts
    data = {
        'Motif': motif_labels,
        'Count': motif_counts
    }
    df = pd.DataFrame(data)
    
    # Ensure colors match the order of the motifs in the bar chart
    color_list = [motif_colors[int(label)] for label in sorted(motif_count.keys()) if motif_colors is not None and int(label) < len(motif_colors)]    
    # Plot using Altair, ensuring the colors match the motifs
    bar_chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Motif', sort=None),
        y='Count',
        tooltip=['Motif', 'Count'],
        color=alt.Color('Motif', scale=alt.Scale(domain=motif_labels, range=color_list))  # Ensure the scale matches correctly
    ).properties(
        width=200,
        height=200,
        title="Motif Occurrences"
    )

    bar_chart = bar_chart.configure_axisX(labelAngle=0)
    st.altair_chart(bar_chart, use_container_width=True)



def create_motif_dataframe(sequences, motif_colors, motif_ids, spans_list, motif_names):
    data = []
    interruptions_dict = set()
    for idx, sequence in enumerate(sequences):
        sequence_name = sequence['name']
        motif_ids_seq = motif_ids[idx]
        spans = spans_list[idx]
        ranges = parse_motif_range(spans)
        sequence_length = len(sequence['sequence'])
        previous_end = 0
        interruptions_dict_sample = {}
        for i, (start, end) in enumerate(ranges):
            motif = motif_ids_seq[i]
            color = motif_colors[int(motif)]

            # Handle interruption between motifs
            if start > previous_end:
                data.append({
                    'Sample': sequence_name,
                    'Start': previous_end,
                    'End': start,
                    'Motif': 'Interruption',
                    'Color': '#FF0000',
                    'Sequence': sequence['sequence'][previous_end:start],
                })
                if sequence['sequence'][previous_end:start] in interruptions_dict_sample:
                    interruptions_dict_sample[sequence['sequence'][previous_end:start]] += 1
                else:
                    interruptions_dict_sample[sequence['sequence'][previous_end:start]] = 1
            # Add motif data
            data.append({
                'Sample': sequence_name,
                'Start': start,
                'End': end + 1,  # Altair works with non-inclusive end
                'Motif': motif_names[int(motif)],
                'Color': color,
                'Sequence': sequence['sequence'][start:end+1],
            })

            previous_end = end + 1
        def len_inturruption_is_equal_to_motif_length(motif_names,k):
            for motif in motif_names:
                if len(k) == len(motif):
                    return True
            return False
        # filter the interruption_dict and keep only the one with sequence == any of the motif_lengths and are more than 1
        inturruptions_dict_sample = {k: v for k, v in interruptions_dict_sample.items() if  len_inturruption_is_equal_to_motif_length(motif_names,k) and v > 1}
        # add the inturruption_dict_sample to the inturruption_dict
        for k,v in inturruptions_dict_sample.items():
            interruptions_dict.add(k)
            

        # Add interruption after the last motif if any
        if previous_end < sequence_length:
            data.append({
                'Sample': sequence_name,
                'Start': previous_end,
                'End': sequence_length,
                'Motif': 'Interruption',
                'Color': '#FF0000',
                'Sequence': sequence['sequence'][previous_end:],
            })
    
    # print the interruptions as inturrption seen : 
    if interruptions_dict:
       
        interruptions_list = " | ".join([f"`{seq}`" for seq in interruptions_dict])
        st.markdown(f"**Interruptions Observed:** {interruptions_list}")
    else:
        st.subheader("No Significant Interruptions Detected")
    return pd.DataFrame(data)


def plot_Cohort_results(cohort_records):
    sequences = []
    span_list = []
    motif_ids_list = []
    # make space between the last print 
    sort_by = st.radio("Sort by:", ("Sample Name" ,"Value" ), horizontal=True, key="sort_by_cohort")
    for key in cohort_records.keys():
        sequences.append({'name': f'{key}_alle1', 'sequence': cohort_records[key]['alt_allele1']})
        span_list.append(cohort_records[key]['spans'][1])
        motif_ids_list.append(cohort_records[key]['motif_ids_h1'])
        if cohort_records[key]['alt_allele2'] != '':
            sequences.append({'name': f'{key}_alle2', 'sequence': cohort_records[key]['alt_allele2']})
            span_list.append(cohort_records[key]['spans'][2])
            motif_ids_list.append(cohort_records[key]['motif_ids_h2'])
        else :
            sequences.append({'name': f'{key}_alle2', 'sequence': ""})
            span_list.append("")
            motif_ids_list.append([0])

    motif_names = cohort_records[list(cohort_records.keys())[0]]['motifs']
    record = cohort_records[list(cohort_records.keys())[0]]
    motif_colors, df = stack_plot(record, motif_names, sequences, span_list, motif_ids_list,sort_by)

    # Filterung der Daten
    figure = go.Figure()

    # Einzigartige Proben erhalten
    unique_samples = df['Sample'].unique()
    # Unterbrechungen entfernen
    unique_samples = [sample for sample in unique_samples if sample != "Interruption"]
    
    # Scatter-Plots für jede Probe und jedes Motiv hinzufügen
    ref_data = []
    allele_data = []

    # Scatter-Plots für jede Probe und jedes Motiv hinzufügen
    
    for sample in unique_samples:
        sample_df = df[df['Sample'] == sample]
        unique_motifs = sample_df['Motif'].unique()
        unique_motifs = [motif for motif in unique_motifs if motif != "Interruption"]
        for motif in unique_motifs:
            motif_df = sample_df[sample_df['Motif'] == motif]
            trace = go.Scatter(
                x=[motif],
                y=[len(motif_df)],
                mode='markers',
                name=sample,
                marker=dict(size=20)
            )
            # Daten für "Ref" und "Allel" sammeln
            if sample in ["Ref", "Allel1", "Allel2"]:
                if sample == "Ref":
                    # change the marker to X
                    trace['marker']['symbol'] = "x"
                    ref_data.append(trace)
                else:
                    trace['marker']['symbol'] = "triangle-down"
                    
                    allele_data.append(trace)
            else:
                figure.add_trace(trace)

    # Daten für "Ref" und "Allel" am Ende hinzufügen
    for trace in ref_data + allele_data:
        figure.add_trace(trace)
    import numpy as np
    # Spezifische Farben für Referenz, Allele und HGSVC definieren
    color_mapping = { sample: np.random.choice(px.colors.qualitative.Plotly) for sample in unique_samples}
    color_mapping["Ref"] = "green"

    # Trace-Farben basierend auf dem Probenamen aktualisieren
    figure.for_each_trace(lambda trace: trace.update(marker=dict(color=color_mapping.get(trace.name, 'gray'))))

    # Doppelte Legendeneinträge entfernen und nur die definierten in color_mapping behalten
    unique_legend_names = set()
    for trace in figure.data:
        if trace.name in unique_legend_names:
            trace.showlegend = False
        elif trace.name in color_mapping:
            unique_legend_names.add(trace.name)
            trace.showlegend = True
        else:
            trace.showlegend = False
    # add the colore gray to the legend and call it HGSVC
    
                # Layout mit Titel und Achsenbeschriftungen aktualisieren
    figure.update_layout(
        title="Motif Occurrences",
        xaxis_title="Motif",
        yaxis_title="Count",
    )

    # X-Achsen-Motive basierend auf ihrer Farbe färben
    xaxis_colors = {motif: motif_colors[idx] for idx, motif in enumerate(motif_names)}
    figure.update_xaxes(tickmode='array', tickvals=list(xaxis_colors.keys()), ticktext=[
        f'<span style="color:{xaxis_colors[motif]}">{motif}</span>' for motif in xaxis_colors.keys()
    ], tickangle=45)
    # show the y axis from 0 to the maximum number of motifs
    figure.update_yaxes(range=[0, df['Sample'].value_counts().max()])

    # Plotly-Figur in Streamlit anzeigen
    st.plotly_chart(figure, use_container_width=True)
    bar_plot_motif_count(df, sort_by=sort_by)
    # plot a heatmap of the count of the motifs for each sample
    # Pivot the DataFrame to get the count of motifs for each sample
    motif_counts = df[df['Motif'] != 'Interruption'].groupby(['Sample', 'Motif']).size().reset_index(name='Count')

    # Create a pivot table for the heatmap
    heatmap_data = motif_counts.pivot(index='Sample', columns='Motif', values='Count').fillna(0)

    # Reset index and melt the DataFrame
    heatmap_data_long = heatmap_data.reset_index().melt(id_vars='Sample', var_name='Motif', value_name='Count')

    # Create a heatmap using Altair
    plot_heatmap(heatmap_data_long, sort_by)

    
def bar_plot_motif_count(df, sort_by="Value"):
    df = df[df['Motif'] != "Interruption"]
    
    # Calculate total copy number for each sample across all motifs
    total_copy_number = df.groupby('Sample').size().reset_index(name='Total Copy Number')
    total_copy_number = total_copy_number.sort_values(by='Total Copy Number', ascending=False)

    # Sort by either value or name
    if sort_by == "Value":
        x_sort = alt.EncodingSortField(field='Total Copy Number', op='sum', order='descending')
    else:
        x_sort = alt.SortField(field='Sample', order='ascending')



    unique_samples = list(total_copy_number['Sample'].apply(lambda x: x.rsplit('_', 1)[0]).unique())
    color_palette =   px.colors.qualitative.Vivid + px.colors.qualitative.Safe + px.colors.qualitative.Dark24 + px.colors.qualitative.Prism 
    

    color_mapping = {sample: color_palette[i] for i, sample in enumerate(unique_samples)}
    color_mapping = {sample: color_mapping[sample.rsplit('_', 1)[0]] for sample in total_copy_number['Sample']}
    # make the colors of both alleles the same

    # Create the bar chart
    bar_chart = alt.Chart(total_copy_number).mark_bar().encode(
        x=alt.X('Sample', sort=x_sort),
        y='Total Copy Number',
        tooltip=['Sample', 'Total Copy Number'],
        # Use Sample for consistent coloring across alleles
        color=alt.Color('Sample', scale=alt.Scale(domain=list(color_mapping.keys()), range=list(color_mapping.values())))
    ).properties(
        width=600,
        height=400,
        title="Total Copy Number per Sample"
    )
    
    # Remove the legend
    bar_chart = bar_chart.configure_legend(orient='none', disable=True)

    # Display bar chart in Streamlit
    st.altair_chart(bar_chart, use_container_width=True)

    

def plot_HGSVC_VS_allele(record, hgsvc_records, motif_names):
    sequences = []
    span_list = []
    motif_ids_list = []
    sort_by = st.radio("Sort by:", ("Value", "Sample Name"), horizontal=True, key="sort_by_HGSVC")
    sequences.append({'name': "Ref", 'sequence': record['ref_allele']})
    span_list.append(record['spans'][0])
    motif_ids_list.append(record['motif_ids_ref'])
    sequences.append({'name': "Allel1", 'sequence':  record['alt_allele1']})
    span_list.append(record['spans'][1])
    motif_ids_list.append(record['motif_ids_h1'])
    if record['alt_allele2'] != '' and record['spans'][2] != None:
        sequences.append({'name': "Allel2", 'sequence': record['alt_allele2']})
        span_list.append(record['spans'][2])
        motif_ids_list.append(record['motif_ids_h2'])

    for key in hgsvc_records.keys():
        if hgsvc_records[key]['alt_allele'] == '':
            continue
        sequences.append({'name': key, 'sequence': hgsvc_records[key]['alt_allele']})
        span_list.append(hgsvc_records[key]['spans'][1])
        motif_ids_list.append(hgsvc_records[key]['motif_ids_h'])

    motif_colors, df = stack_plot(record, motif_names, sequences, span_list, motif_ids_list,sort_by)
    bar_plot_motif_count(df, sort_by)


    # Filterung der Daten
    figure = go.Figure()

    # Einzigartige Proben erhalten
    unique_samples = df['Sample'].unique()
    # Unterbrechungen entfernen
    unique_samples = [sample for sample in unique_samples if sample != "Interruption"]

    # Scatter-Plots für jede Probe und jedes Motiv hinzufügen
    ref_data = []
    allele_data = []

    # Scatter-Plots für jede Probe und jedes Motiv hinzufügen
    for sample in unique_samples:
        sample_df = df[df['Sample'] == sample]
        unique_motifs = sample_df['Motif'].unique()
        unique_motifs = [motif for motif in unique_motifs if motif != "Interruption"]
        for motif in unique_motifs:
            motif_df = sample_df[sample_df['Motif'] == motif]
            trace = go.Scatter(
                x=[motif],
                y=[len(motif_df)],
                mode='markers',
                name=sample,
                marker=dict(size=20)
            )
            # Daten für "Ref" und "Allel" sammeln
            if sample in ["Ref", "Allel1", "Allel2"]:
                if sample == "Ref":
                    # change the marker to X
                    trace['marker']['symbol'] = "x"
                    ref_data.append(trace)
                else:
                    trace['marker']['symbol'] = "triangle-down"
                    
                    allele_data.append(trace)
            else:
                figure.add_trace(trace)

    # Daten für "Ref" und "Allel" am Ende hinzufügen
    for trace in ref_data + allele_data:
        figure.add_trace(trace)
    # Spezifische Farben für Referenz, Allele und HGSVC definieren
    color_mapping = {
        "Ref": "green",
        "Allel1": "red",
        "Allel2": "orange",
        "HGSVC": "blue"
    }

    # Trace-Farben basierend auf dem Probenamen aktualisieren
    figure.for_each_trace(lambda trace: trace.update(marker=dict(color=color_mapping.get(trace.name, 'gray'))))

    # Doppelte Legendeneinträge entfernen und nur die definierten in color_mapping behalten
    unique_legend_names = set()
    for trace in figure.data:
        if trace.name in unique_legend_names:
            trace.showlegend = False
        elif trace.name in color_mapping:
            unique_legend_names.add(trace.name)
            trace.showlegend = True
        else:
            trace.showlegend = False
    # add the colore gray to the legend and call it HGSVC
    name = "HGSVC" if st.session_state.analysis_mode == "indivisual sample" else st.session_state.analysis_mode 
    figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color='gray', size=20), name=name))
    # Layout mit Titel und Achsenbeschriftungen aktualisieren
    figure.update_layout(
        title="Motif Occurrences",
        xaxis_title="Motif",
        yaxis_title="HGSVC vs Alleles",
    )

    # X-Achsen-Motive basierend auf ihrer Farbe färben
    xaxis_colors = {motif: motif_colors[idx] for idx, motif in enumerate(motif_names)}
    figure.update_xaxes(tickmode='array', tickvals=list(xaxis_colors.keys()), ticktext=[
        f'<span style="color:{xaxis_colors[motif]}">{motif}</span>' for motif in xaxis_colors.keys()
    ], tickangle=45)
    # show the y axis from 0 to the maximum number of motifs
    figure.update_yaxes(range=[0, df['Sample'].value_counts().max()])
    
    # Plotly-Figur in Streamlit anzeigen
    st.plotly_chart(figure, use_container_width=True)
    #Motif Overlap Radial Chart
    # Pivot tables for HGSVC and sample data
    pivot_hgsvc = pd.pivot_table(df[df['Sample'] == 'HGSVC'], index='Motif', columns='Sample', values='Length', aggfunc='count', fill_value=0)
    pivot_sample = pd.pivot_table(df[df['Sample'] != 'HGSVC'], index='Motif', columns='Sample', values='Length', aggfunc='count', fill_value=0)

    # Convert pivot tables to long format for Altair
    pivot_hgsvc_long = pivot_hgsvc.reset_index().melt(id_vars='Motif', var_name='Sample', value_name='Count')
    pivot_sample_long = pivot_sample.reset_index().melt(id_vars='Motif', var_name='Sample', value_name='Count')

    # Combine the data for comparison
    combined_data = pd.concat([pivot_hgsvc_long, pivot_sample_long])

    # Create Altair heatmap
    plot_heatmap(combined_data, sort_by=sort_by)

def plot_heatmap(combined_data, sort_by="Value"):
    

    if sort_by == "Value":
        x_sort = alt.EncodingSortField(field='Count', op='sum', order='descending')
        y_sort = alt.EncodingSortField(field='Count', op='sum', order='descending')
    else:
        x_sort = alt.SortField(field='Sample', order='ascending')
        y_sort = alt.SortField(field='Motif', order='ascending')

    heatmap = alt.Chart(combined_data).mark_rect().encode(
        x=alt.X('Sample:N', title='Sample', sort=x_sort),
        y=alt.Y('Motif:N', title='Motif', sort=y_sort),
        color=alt.Color('Count:Q', scale=alt.Scale(scheme='reds'), title='Count'),
        tooltip=['Sample', 'Motif', 'Count']
    ).properties(
        width=600,
        height=400,
        title='Motif Occurrences Heatmap'
    )
    # make the legend above the plot 
    heatmap = heatmap.configure_legend(orient='top')
    # Display heatmap in Streamlit
    st.altair_chart(heatmap, use_container_width=True)

def stack_plot(record, motif_names, sequences, span_list, motif_ids_list, sort_by="Value"):
    motif_colors = get_color_palette(len(record['motifs']))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
    
    region  = record['chr'] + ":" + str(record['pos']-1) + "-" + str(record['stop']-1)

    df = create_motif_dataframe(sequences, motif_colors, motif_ids_list, span_list, motif_names)

    df['Length'] = df['End'] - df['Start']
    
    default_hight = 600 
    chart_height = max(default_hight, len(sequences) * 7)

    df['Sample'] = df['Sample'].apply(lambda x: x.replace("_pathogenic", ""))
    df['Sample'] = pd.Categorical(df['Sample'], categories=sorted(df['Sample'].unique()), ordered=True)
    df = df.sort_values(by=['Sample', 'Start', 'End'])
    df['Order'] = range(len(df))

    df_filtered = df[df['Motif'] != 'Interruption']

    for sample in df_filtered['Sample'].unique():
        df.loc[df['Sample'] == sample, 'Total Copy Number'] = df[df['Sample'] == sample].shape[0]
    
    min_copy_number = df['Total Copy Number'].min()
    max_copy_number = df['Total Copy Number'].max()

    pathogenic_thresold = 0
    if region in st.session_state.pathogenic_TRs:
        
        pathogenic_thresold = st.session_state.pathogenic_TRs[region]['pathogenicity_threshold']
        pathogenic_thresold = min(pathogenic_thresold) if pathogenic_thresold != [] else 0
        gene_name = st.session_state.pathogenic_TRs[region]['gene']
        st.markdown(f"""<div style="height: 20px;"></div>""", unsafe_allow_html=True)
        # Genname anzeigen
        st.markdown(f"""
            <div style="display: flex; justify-content: center; font-size: 20px; color: #4CAF50;">
            <strong>Gene Name:</strong> {gene_name}
            </div>
        """, unsafe_allow_html=True)
        # show the pathogenicity threshold 
        if pathogenic_thresold == 0:
            pathogenic_thresold_print = "Unknown"
        else:
            pathogenic_thresold_print = pathogenic_thresold
        st.markdown(f"""
            <div style="display: flex; justify-content: center; font-size: 20px; color: #FF5733;">
                <strong>Pathogenicity Threshold: </strong>  {pathogenic_thresold_print}
            </div>
        """, unsafe_allow_html=True)
    st.markdown(f"""
        <div style="display: flex; justify-content: space-between; font-size: 20px; color: #FF5733;">
            <div>
                <strong>Minimum Copy Number:</strong> {min_copy_number}
            </div>
            <div>
                <strong>Maximum Copy Number:</strong> {max_copy_number}
            </div>
        </div>
    """, unsafe_allow_html=True)

    if sort_by == "Value":
        y_sort = alt.EncodingSortField(field='Length', op='sum', order='descending')
    else:
        y_sort = alt.SortField(field='Sample', order='ascending')
    
    df["Total Copy Number"].fillna(0, inplace=True)

    df["pathogenic"] = df.apply(lambda row: "Pathogenic" if  (row["Motif"] != "Interruption" and row["Total Copy Number"] >= pathogenic_thresold  ) else "Not Pathogenic", axis=1)
    chart = alt.Chart(df).mark_bar().encode(
        y=alt.Y(
            'Sample', 
            sort=y_sort, 
            axis=alt.Axis(labelOverlap=False, ticks=False)  
        ),
        x=alt.X('Length', title='Length', stack='zero'),
        color=alt.Color('Motif', scale=alt.Scale(domain=list(motif_names) + ['Interruption'], range=list(motif_colors.values()) + ['#FF0000'])),
        order=alt.Order('Order', sort='ascending'),
        tooltip=['Sample', 'Motif', 'Start', 'End', 'Sequence','pathogenic']
    ).properties(
        width=800,
        height=chart_height,
        title="Motif Occurrences"
    ).interactive()

    if chart_height > default_hight:
        chart = chart.configure_axisY(labelFontSize=0)
    # check if the pathogenicity threshold is greater than 0 and if all the samples are not pathogenic
    pathogenic_thresold_length = pathogenic_thresold * len(motif_names[0])
    # Abstand hinzufügen

        
    if pathogenic_thresold > 0:

        # Pathogenicity threshold multiplizieren
        
        if  df.groupby('Sample')['Length'].sum().max() > pathogenic_thresold_length:
            # Regel hinzufügen und beschriften
            rule = alt.Chart(pd.DataFrame({'x': [pathogenic_thresold_length], 'label': ['Pathogenic Threshold']})).mark_rule(color='red').encode(
                x='x',
                tooltip=['label', 'x'],
            )
            
            # make the line thicker
            rule = rule.encode(size=alt.value(5))
            
            # Diagramme kombinieren
            combined_chart = alt.layer(chart, rule).resolve_scale(y='independent')
            

            # Konfiguration anwenden
            combined_chart = combined_chart.configure_axis(
                labelFontSize=10,
                titleFontSize=12,
            )

            combined_chart = combined_chart.properties(
                padding={'left': 10, 'right': 50, 'top': 30, 'bottom': 10}
            )

            # for all the samples that cross the pathogenic threshold color their name on the y-axis red
        

            st.altair_chart(combined_chart, use_container_width=True)
        else:
            st.altair_chart(chart, use_container_width=True)
    else:
        st.altair_chart(chart, use_container_width=True)
    return motif_colors, df


# Function to visualize tandem repeat with highlighted motifs on the sequence
def visulize_TR_with_dynamic_sequence(record,hgsvc_records, left_column, right_column,motif_colors,CN1_col,CN2_col, show_comparison):
    motif_names = record['motifs']
    reference_copy_number = record['ref_CN']
    motif_count_ref = count_motifs(record['motif_ids_ref'], record['spans'][0])
    found_motifs_ref = list(motif_count_ref.keys())
    found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]

    motif_count_h1 = count_motifs(record['motif_ids_h1'], record['spans'][1])
    found_motifs_h1 = list(motif_count_h1.keys())
    found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
    motif_count_h2 = {}
    motif_count_h2 = {}
    found_motifs_h2 = []
    if record['alt_allele2'] != '' and record['spans'][2] != None:
        motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][2])
        found_motifs_h2 = list(motif_count_h2.keys())
        found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
        motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
    total_copy_number_h1 = str(record['spans'][1]).count('-')
    CN1_col.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 1 Total copy number:</strong> {total_copy_number_h1}
            </div>
        """, unsafe_allow_html=True)
    if record['alt_allele2'] != '' and record['spans'][2] != None:
            total_copy_number_h2 = str(record['spans'][2]).count('-')
            CN2_col.markdown(f"""
                <div style="font-size: 20px; color: #FF5733;">
                    <strong>Allele 2 Total copy number:</strong> {total_copy_number_h2}
                </div>
            """, unsafe_allow_html=True) 
    else: 
        CN2_col.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 2 Total copy number:</strong> &lt;DEL&gt;
            </div>
        """, unsafe_allow_html=True)
    with right_column:
        display_motif_legend(motif_names, motif_colors, right_column)
    with left_column:
        tab1, tab2,tab3 = st.tabs(["Alleles", "Alleles vs Ref", "Alleles vs Pop"])
        # create a spot for the plots to refresh
        with tab2:
            display_dynamic_sequence_with_highlighted_motifs("Ref", record['ref_allele'], record['motif_ids_ref'], record['spans'][0], motif_colors, motif_names)
            alt_allele1 = record['alt_allele1']
            if record['alt_allele1'] != '.' and record['spans'][1] != None:
                display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names)
            if record['alt_allele2'] != '.' and record['spans'][2] != None:
                alt_allele2 = record['alt_allele2'] #if record['alt_allele2'] != "." else record['ref_allele']
                display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names)
        with tab1:
        # Render the scrollable sequence with highlighted motifs for allele 1
            alt_allele1 = record['alt_allele1']
            if record['alt_allele1'] != '.' and record['spans'][1] != None:
                display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names)

            # Create an empty container for the plot to refresh
            plot_container_h1 = st.empty()
            with plot_container_h1:
                if show_comparison == False:
                    plot_motif_bar(motif_count_h1, motif_names, motif_colors)

            if record['alt_allele2'] != '' and record['spans'][2] != None:
                alt_allele2 = record['alt_allele2'] #if record['alt_allele2'] != "." else record['ref_allele']
                display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names)
                
                # Create another empty container for the plot to refresh
                plot_container_h2 = st.empty()
                with plot_container_h2:
                    if show_comparison == False:
                        plot_motif_bar(motif_count_h2, motif_names, motif_colors)
        with tab3:
            # get the records for the population
            # add the alleles 
            plot_HGSVC_VS_allele(record, hgsvc_records, motif_names)


def display_summary_statistics(records):
    total_motifs = 0


    for record in records.values():
        total_motifs += len(record['motifs'])
      

    st.sidebar.markdown("### Summary Statistics")
    st.sidebar.text(f"Total number Records: {len(records)}")
    st.sidebar.text(f"Total Motifs: {total_motifs}")

def count_motifs(motif_ids, spans):
    ranges = parse_motif_range(spans)
    motif_count = {}
    
    for idx, motif in enumerate(motif_ids):
        if motif in motif_count:
            motif_count[motif] += 1
        else:
            motif_count[motif] = 1
    
    return motif_count


def display_motif_legend(motifs, motif_colors, right_column):
    st.markdown("### Motif Legend")
    st.markdown('<div style="max-height:400px; overflow-y:scroll;">', unsafe_allow_html=True)  # Scrollable container
    if  isinstance(motifs, tuple):
        motifs = list(motifs)
    elif not isinstance(motifs, list):
        motifs = [motifs]
    for idx, motif in enumerate(motifs):
        color = motif_colors[idx]
        # Assign each legend item a unique class for hover effect
        st.markdown(
            f'<div id="legend-motif-{idx}" class="legend-item motif-{idx}" style="background-color:{color};color:white;padding:5px;margin-bottom:10px;border-radius:5px;'
            f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);'
            f' white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="{motif}">'
            f' Motif {idx}: {motif}</div>', unsafe_allow_html=True)
    # show the inturruptions as well as the gray catagory
    st.markdown(
        f'<div id="legend-motif-interruption" class="legend-item motif-interruption" style="background-color:#FF0000;color:black;padding:5px;margin-bottom:10px;border-radius:5px;'
        f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
        f' Interruption</div>', unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)

@st.cache_data()
def get_records_info(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
         
    idx = 0 
    cohorts_map = {}
    for rec in vcf:
        cohorts_map[idx] = rec.id
        idx += 1
    return cohorts_map
            
def fetch_vcf_region(vcf_file_path, region):
    vcf = pysam.VariantFile(vcf_file_path)
    records = []

    for rec in vcf.fetch(region=region):  # Fetch based on user input region
        ids_h1 = rec.info['MOTIF_IDs_H1']
        ids_h2 = rec.info['MOTIF_IDs_H2']
        alt_allele1 = "."
        alt_allele2 = "."
        ref_allele = rec.alts[0] if rec.alts != None else ""
        print(rec.alts)
        if rec.alts != None:
            if len(rec.alts) > 0:
                alt_allele1 = rec.alts[0]
                if alt_allele1 == '.':
                    alt_allele1 = ''
            if len(rec.alts) > 1:
                alt_allele2 = rec.alts[1]
            else:
                if ids_h1 == ids_h2:
                    alt_allele2 = alt_allele1
                alt_allele2 = ''
        if ids_h2 is None:
            ids_h2 = []
        motif_names = rec.info['MOTIFS']
        if isinstance(motif_names, tuple):
            motif_names = list(motif_names)
        elif not isinstance(motif_names, list):
            motif_names = [motif_names]

        record = {
            'chr': rec.chrom,
            'pos': rec.pos,
            'motifs': motif_names,
            'motif_ids_h1': ids_h1,
            'motif_ids_h2': ids_h2,
            'spans': rec.samples[0]['SP'],
            'ref_allele': ref_allele,
            'alt_allele1': alt_allele1,
            'alt_allele2': alt_allele2
        }
        records[rec.id] = record
    return records


# Streamlit UI
st.sidebar.markdown("""
    <style>
        :root {
            --primary-color: #4CAF50;
            --secondary-color: #333;
            --background-color: #f9f9f9;
            --text-color: #555;
            --border-color: #4CAF50;
            --hover-color: #388E3C;
        }

        @media (prefers-color-scheme: dark) {
            :root {
                --primary-color: #90EE90;
                --secondary-color: #ddd;
                --background-color: #333;
                --text-color: #ccc;
                --border-color: #90EE90;
                --hover-color: #76C76A;
            }
        }

        .sidebar-container {
            text-align: center;
            font-family: 'Segoe UI', Tahoma, sans-serif;
        }

        .sidebar-container h1 {
            color: var(--primary-color);
            font-size: 28px;
            margin-bottom: 10px;
            animation: fadeIn 1.5s ease-in-out;
        }

        .sidebar-container hr {
            border: 1px solid var(--border-color);
            margin: 20px 0;
            animation: grow 1s ease-in-out;
        }

        .sidebar-container p {
            color: var(--text-color);
            font-size: 16px;
            margin-bottom: 20px;
            animation: fadeIn 2s ease-in-out;
        }

        .upload-container, .instructions-container {
            background-color: var(--background-color);
            padding: 5px;
            border-radius: 30px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            margin-bottom: 20px;
            transition: transform 0.3s ease;
        }

        .upload-container:hover, .instructions-container:hover {
            transform: translateY(-5px);
        }

        .upload-container h3, .instructions-container h3 {
            color: var(--primary-color);
            text-align: center;
            margin-bottom: 10px;
        }

        .upload-container p, .instructions-container ul {
            color: var(--secondary-color);
            font-size: 14px;
            text-align: center;
        }

        .button-container {
            text-align: center;
        }

        .button-container button {
            background-color: var(--primary-color);
            color: white;
            padding: 10px 20px;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            transition: background-color 0.3s;
        }

        .button-container button:hover {
            background-color: var(--hover-color);
        }

        @keyframes fadeIn {
            0% { opacity: 0; }
            100% { opacity: 1; }
        }

        @keyframes grow {
            0% { width: 0; }
            100% { width: 100%; }
        }
    
    </style>
""", unsafe_allow_html=True)

# Sidebar content
st.sidebar.markdown("""
    <style>

        @media (prefers-color-scheme: light) {
            :root {
                --sidebar-text-color: darkgray;
            }
        }
    </style>
    <div class='sidebar-container' style="color: var(--sidebar-text-color);">
        <h1 style="color: var(--sidebar-text-color);">Tandem Repeat Visualization</h1>
        <hr>
        <p style="color: var(--sidebar-text-color);">Visualize and analyze tandem repeats in your genomic data.</p>
    </div>
""", unsafe_allow_html=True)

# Container for file upload make it smaller
# st.sidebar.markdown("""
#     <div class='upload-container'>
#         <h3>Upload VCF File</h3>
#         <p>Please enter the path of your VCF file to get started.</p>
        
#     </div>
# """, unsafe_allow_html=True)




path_changed = False
old_vcf_file_path = st.session_state.get('vcf_file_path', None)


# Check if the analysis mode has changed
previous_analysis_mode = st.session_state.get('previous_analysis_mode', None)
st.session_state.analysis_mode = st.sidebar.radio("Select the type of analysis", ("indivisual sample", "Cohort"))

# If the analysis mode has changed, reset the index to 0
if previous_analysis_mode != st.session_state.analysis_mode:
    st.session_state.regions_idx = 0

# Store the current analysis mode for future comparison
st.session_state.previous_analysis_mode = st.session_state.analysis_mode

if st.session_state.analysis_mode == "indivisual sample":
    st.sidebar.text_input("Enter the path of your VCF file", key="vcf_file_path",value= None, help="Enter the path of your VCF file, the file should be zipped and indexed with tabix")
    if st.session_state.get('vcf_file_path', None) is None:
        st.stop()
    if st.session_state.get('vcf_file_path', None) != old_vcf_file_path:
        path_changed = True
        st.session_state.vcf_file_path = st.session_state.get('vcf_file_path', None)
        st.session_state.pop('records', None)
        st.session_state.pop('records_map', None)
        

    if 'records' not in st.session_state or path_changed:
        if st.session_state.analysis_mode == "indivisual sample":
        
        # posistion the button in the center
            _, _,middle, _ = st.sidebar.columns([1,0.3, 2, 1])
            
            st.session_state.vcf_status = st.sidebar.radio("Select the type of VCF file", ("Healthy", "Pathogenic"))
        if middle.button("Upload VCF File"):
            # try:
            if 'records' not in st.session_state:
                st.session_state.records,st.session_state.records_map = parse_vcf( st.session_state.get('vcf_file_path', None))
                st.session_state.hgsvc_path = "/confidential/tGenVar/Lion/tandemtwister-vis/HGSVC_vcf/"
                if st.session_state.vcf_status == "Pathogenic":
                    st.session_state.file_paths = [f for f in os.listdir(st.session_state.hgsvc_path) if f.endswith('pathogenic.vcf.gz')]
                elif st.session_state.vcf_status == "Healthy":
                    st.session_state.file_paths = [f for f in os.listdir(st.session_state.hgsvc_path) if f.endswith('.vcf.gz')]
                st.session_state.files = [load_vcf(st.session_state.hgsvc_path + f) for f in st.session_state.file_paths]
                
                

    subheader = st.empty()

    if 'records_map' in st.session_state:
        if 'regions_idx' not in st.session_state:
            st.session_state.regions_idx = 0

        # Sidebar for region navigation
        st.sidebar.markdown("### Select Region to Visualize")
        # activate the variable when pressing inter button
        region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None, key="region", help="Enter the region in the format: chr:start-end")
        #_level = st.sidebar.slider('Zoom Level', min_value=1, max_value=100, value=100)
        
        display_option = st.sidebar.radio("Select Display Type", 
                                    ( "Bars","Sequence with Highlighted Motifs"))

        col1, middel, col2 = st.columns([1.5,3, 1])  # Adjust the ratio [1, 1] to control spacing between buttons
        REF, CN1_col, CN2_col = st.columns([1, 1, 1])
        # Place the "Previous region" and "Next region" buttons in these columns
        with col1:
            if st.button("Previous region"):
                region = None
                st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

        with col2:
            if st.button("Next region"):
                region = None
                st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.records_map) - 1)
   
        if region and region != st.session_state.get('previous_region', None):
            st.session_state.previous_region = region
            try:
                chr_input, start_end_input = region.split(':')
                start_input, end_input = map(int, start_end_input.split('-'))
                # start_input
                # end_input-=

                input_region = f"{chr_input}:{start_input}-{end_input}"
                record_key = st.session_state.records[input_region]
                # update the region index based on the record key
                st.session_state.regions_idx = list(st.session_state.records_map.values()).index(input_region)
            

            except:
                try:
                    chr_input, start_input, end_input = re.split(r'\s+', region)
                    start_input, end_input = int(start_input), int(end_input)
                    input_region = f"{chr_input}:{start_input}-{end_input}"
                    record_key = st.session_state.records[input_region]
                    # update the region index based on the record key
                    st.session_state.regions_idx = list(st.session_state.records_map.values()).index(input_region)
                except:
                    st.sidebar.info("Invalid region format, showing the first record")
                    record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]
        else:
            record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]


        record = parse_record( st.session_state.get('vcf_file_path', None), record_key)
        hgsvc_records = get_results_hgsvc_pop(record_key, st.session_state.files ,st.session_state.file_paths)
        #st.session_state.hgsvc_pop_records, = parse_hgsvc_pop(vcf_status,region)
        if len(record["motif_ids_h1"]) == 0 and len(record["motif_ids_h2"]) == 0:
            st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
            st.stop()
            
        middel.markdown(f"""
                <div id="tandem-repeat-region" class="region-container" style="font-size: 25px; margin-bottom: 10px;">
                    <strong>Tandem Repeat Region: {st.session_state.records_map[st.session_state.regions_idx]}</strong>
                </div>
            """, unsafe_allow_html=True)

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
        # show the reference copy number
        REF.markdown(f"""
            <div style="font-size: 20px; color: #4CAF50; margin-bottom: 10px;">
                <strong>Reference Copy Number:</strong> {record['ref_CN']}
            </div>
        """, unsafe_allow_html=True)

        left_column, right_column = st.columns([4, 1])
        # define the motif colors
        motif_colors = get_color_palette(len(record['motifs']))
        motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
    
        col1,col2 = st.sidebar.columns([1,1])



        if display_option == "Sequence with Highlighted Motifs":
            visulize_TR_with_dynamic_sequence(record,hgsvc_records, left_column, right_column,motif_colors,CN1_col,CN2_col, st.session_state.get('show_comparison', False))

        elif display_option == "Bars":
            display_motifs_with_bars(record, left_column, right_column,motif_colors,CN1_col,CN2_col, st.session_state.get('show_comparison', False))
elif st.session_state.analysis_mode == "Cohort":
    cohort_path = st.sidebar.text_input("Enter the path to the cohort results", value=None, key="cohort_path")
    # if the path doesn't end with a slash add it
    if cohort_path is not None and not cohort_path.endswith('/'):
        cohort_path += '/'
    st.session_state.path_to_cohort = cohort_path
    if st.session_state.path_to_cohort is None:
        st.stop()
    if st.sidebar.button("Load Cohort"):
        st.session_state.cohort_file_paths = [f for f in os.listdir(st.session_state.path_to_cohort) if f.endswith('.vcf.gz')]
        st.session_state.cohort_files = [load_vcf(st.session_state.path_to_cohort + f) for f in st.session_state.cohort_file_paths]
        st.session_state.cohorts_records_map = get_records_info(st.session_state.path_to_cohort+ st.session_state.cohort_file_paths[0])
    col1, middel, col2 = st.columns([1.5,3, 1])  # Adjust the ratio [1, 1] to control spacing between buttons
        # Place the "Previous region" and "Next region" buttons in these columns
    if 'cohorts_records_map' in st.session_state:
        region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None, key="region", help="Enter the region in the format: chr:start-end")
        


        
        if 'regions_idx' not in st.session_state:
            st.session_state.regions_idx = 0
        with col1:
            if st.button("Previous region"):
                region = None
                st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

        with col2:
            if st.button("Next region"):
                region = None
                st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len( st.session_state.cohorts_records_map )-1)
        if region and region != st.session_state.get('previous_region', None):
            st.session_state.previous_region = region
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

            region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
        mode_placeholder = st.empty()

        # JavaScript code to detect dark or light mode and return the value
        # check the background color of the page
        
        # print the region name in the middle 

        middel.markdown(f"""
            <div id="tandem-repeat-region" class="region-container" style="font-size: 25px; margin-bottom: 10px;">
                <strong>Tandem Repeat Region: {region}</strong>
            </div>
        """, unsafe_allow_html=True)

        # Inject the CSS to handle light and dark mode
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

        st.session_state.cohort_results = get_results_cohort(region, st.session_state.cohort_files, st.session_state.cohort_file_paths)
        if 'cohort_results' in st.session_state:   
            region = st.session_state.regions_idx
            
            plot_Cohort_results(st.session_state.cohort_results)
        else:
            st.stop()
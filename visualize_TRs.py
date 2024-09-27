import streamlit as st
import pysam
import re
import tempfile
import random
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib.colors import ListedColormap
import plotly.express as px
import plotly.graph_objects as go
import altair as alt
from collections import OrderedDict

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

def parse_record(vcf_file,region):
    
    vcf = pysam.VariantFile(vcf_file)
    rec = vcf.fetch(region=region)
    # get the record with the id
    for rec in vcf.fetch(region=region):
        break
    ids_h1 = rec.info['MOTIF_IDs_H1']
    ids_h2 = rec.info['MOTIF_IDs_H2']
    alt_allele1 = "."
    alt_allele2 = "."
    ref_allele = rec.alts[0] if rec.alts != None else ""
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
    if ids_h1 is None:
        ids_h1 = []
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
    return record

def parse_motif_range(motif_range):
    pattern = re.compile(r'\((\d+)-(\d+)\)')
    matches = pattern.findall(motif_range)
    ranges = [(int(start)-1, int(end)-1) for start, end in matches]
    return ranges

# Function to generate a list of n visually distinct colors using a color map
def get_color_palette(n):
    cmap = plt.get_cmap('tab20')  # You can try 'viridis', 'plasma', 'Set1', etc.
    colors = [cmap(i) for i in range(n)]
    return ['#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in colors]

def display_dynamic_sequence_with_highlighted_motifs(sequence, motif_ids, spans, motif_colors, motif_names):
   
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
    
    # Add scrollable container for the sequence display
    st.markdown(f"""
        <div style="font-family:monospace; font-size:16px; width:100%; max-height:120px; overflow-x:auto; white-space:nowrap; padding:10px; border:1px solid black; border-radius:8px;">
            {highlighted_sequence}
        </div>
    """, unsafe_allow_html=True)


# Function to display motifs relative to the sequence, with spans aligned horizontally
def display_motifs_with_bars(record, left_column, right_column,motif_colors):
    motif_names = record['motifs']

    motif_count_h1 = count_motifs(record['motif_ids_h1'], record['spans'][0])
    found_motifs_h1 = list(motif_count_h1.keys())
    found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
    motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][1])
    found_motifs_h2 = list(motif_count_h2.keys())
    found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
    motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
    # iterate over the motifs and set them to 0 if they are not in the motif_count

    # for motif in motif_names:
    #     motif_index = motif_names.index(motif)
    #     if motif not in found_motifs_h1:
    #         motif_count_h1[motif_index] = 0
    #     if motif not in found_motifs_h2:
    #         motif_count_h2[motif_index] = 0

    with right_column:
        display_motif_legend(motif_names, motif_colors, right_column)

    with left_column:

        display_motifs_as_bars(motif_colors, record['motif_ids_h1'], record['spans'][0], record['alt_allele1'], motif_names)
        plot_container_h1 = st.empty()
        with plot_container_h1:
            plot_motif_bar(motif_count_h1, motif_names, motif_colors)
        
        if record['alt_allele2'] != '':
            display_motifs_as_bars(motif_colors, record['motif_ids_h2'], record['spans'][1], record['alt_allele2'], motif_names)
            plot_container_h2 = st.empty()
            with plot_container_h2:
                plot_motif_bar(motif_count_h2, motif_names, motif_colors)


def display_motifs_as_bars(motif_colors, motif_ids, spans, sequence, motif_names):
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

    # Render the motif bars in Streamlit
    st.markdown(bar_container, unsafe_allow_html=True)

def plot_motif_bar(motif_count, motif_names, motif_colors=None):
    # Convert motif indices to names and store counts
    motif_labels = []
    motif_counts = []
    for label, value in sorted(motif_count.items()):
        motif_labels.append(motif_names[int(label)])  # Ensure correct mapping of label to name
        motif_counts.append(value)
    
    # Create DataFrame with motif labels and counts
    data = {
        'Motif': motif_labels,
        'Count': motif_counts
    }
    df = pd.DataFrame(data)

    # Ensure colors match the order of the motifs in the bar chart
    color_list = [motif_colors[int(label)] for label in sorted(motif_count.keys())]

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


# Function to visualize tandem repeat with highlighted motifs on the sequence
def visulize_TR_with_dynamic_sequence(record, left_column, right_column,motif_colors):
    motif_names = record['motifs']


    motif_count_h1 = count_motifs(record['motif_ids_h1'], record['spans'][0])
    found_motifs_h1 = list(motif_count_h1.keys())
    found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
    motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][1])
    found_motifs_h2 = list(motif_count_h2.keys())
    found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
    motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
   
    with right_column:
        display_motif_legend(motif_names, motif_colors, right_column)
    with left_column:
        # create a spot for the plots to refresh
        total_copy_number_h1 = str(record['spans'][0]).count('-')
        st.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 1 Total copy number:</strong> {total_copy_number_h1}
            </div>
        """, unsafe_allow_html=True)
        
        # Render the scrollable sequence with highlighted motifs for allele 1
        alt_allele1 = record['alt_allele1'] if record['alt_allele1'] != "." else record['ref_allele']
        display_dynamic_sequence_with_highlighted_motifs(alt_allele1, record['motif_ids_h1'], record['spans'][0], motif_colors, motif_names)
        
        # Create an empty container for the plot to refresh
        plot_container_h1 = st.empty()
        with plot_container_h1:
            
            plot_motif_bar(motif_count_h1, motif_names, motif_colors)

        if record['alt_allele2'] != '':
            total_copy_number_h2 = str(record['spans'][1]).count('-')
            st.markdown(f"""
                <div style="font-size: 20px; color: #FF5733;">
                    <strong>Allele 2 Total copy number:</strong> {total_copy_number_h2}
                </div>
            """, unsafe_allow_html=True)
            alt_allele2 = record['alt_allele2'] #if record['alt_allele2'] != "." else record['ref_allele']
            display_dynamic_sequence_with_highlighted_motifs(alt_allele2, record['motif_ids_h2'], record['spans'][1], motif_colors, motif_names)
            
            # Create another empty container for the plot to refresh
            plot_container_h2 = st.empty()
            with plot_container_h2:
                plot_motif_bar(motif_count_h2, motif_names, motif_colors)



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


def fetch_vcf_region(vcf_file_path, region):
    vcf = pysam.VariantFile(vcf_file_path)
    records = []
    for rec in vcf.fetch(region=region):  # Fetch based on user input region
        ids_h1 = rec.info['MOTIF_IDs_H1']
        ids_h2 = rec.info['MOTIF_IDs_H2']
        alt_allele1 = "."
        alt_allele2 = "."
        ref_allele = rec.alts[0] if rec.alts != None else ""
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

def lazy_load_vcf(vcf_file_path, start_index, chunk_size):
    vcf = pysam.VariantFile(vcf_file_path)
    records = []
    for i, rec in enumerate(vcf.fetch()):
        if start_index <= i < start_index + chunk_size:
            records.append(rec)
        if i >= start_index + chunk_size:
            break
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
        }

        @media (prefers-color-scheme: dark) {
            :root {
                --primary-color: #90EE90;
                --secondary-color: #ddd;
                --background-color: #333;
                --text-color: #ccc;
                --border-color: #90EE90;
            }
        }

        .sidebar-container {
            text-align: center;
            font-family: Arial, sans-serif;
        }

        .sidebar-container h1 {
            color: var(--primary-color);
        }

        .sidebar-container hr {
            border: 1px solid var(--border-color);
        }

        .sidebar-container p {
            color: var(--text-color);
            font-size: 16px;
        }

        .upload-container {
            background-color: var(--background-color);
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            margin-bottom: 20px;
        }

        .upload-container h3 {
            color: var(--primary-color);
            text-align: center;
        }

        .upload-container p {
            color: var(--secondary-color);
            font-size: 14px;
            text-align: center;
        }

        .instructions-container {
            background-color: var(--background-color);
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        }

        .instructions-container h3 {
            color: var(--primary-color);
            text-align: center;
        }

        .instructions-container ul {
            color: var(--secondary-color);
            font-size: 14px;
        }
    </style>
""", unsafe_allow_html=True)

# Sidebar content
st.sidebar.markdown("""
    <div class='sidebar-container'>
        <h1>Tandem Repeat Visualization</h1>
        <hr>
        <p>Visualize and analyze tandem repeats in your genomic data.</p>
    </div>
""", unsafe_allow_html=True)

# Container for file upload
st.sidebar.markdown("""
    <div class='upload-container'>
        <h3>Upload VCF File</h3>
        <p>Please enter the path of your VCF file to get started.</p>
    </div>
""", unsafe_allow_html=True)
        
# enter the path to the vcf file
path_changed = False
old_vcf_file_path = st.session_state.get('vcf_file_path', None)
# make a file uploader
st.sidebar.text_input("Enter the path to the VCF file", value=None, key="vcf_file_path")
vcf_file_path = st.session_state.get('vcf_file_path', None)

if vcf_file_path != old_vcf_file_path:
    path_changed = True
    st.session_state.vcf_file_path = vcf_file_path
    st.session_state.pop('records', None)
    st.session_state.pop('records_map', None)
# check if records are in st.session_state and if path has changed

if 'records' not in st.session_state or path_changed:
    if st.sidebar.button("upload_file"):
        st.write("uploading file")
        #vcf_file_path = "/confidential/tGenVar/Lion/TandemTwist/MC/HG002_tandemtwister.vcf.gz"
        #vcf_file = st.sidebar.file_uploader("Upload VCF file", type=["vcf", "gz"])
        if 'records' not in st.session_state:
            st.session_state.records,st.session_state.records_map = parse_vcf(vcf_file_path)
      


            
        
        # definde a container for the subheader
    # Example usage in Streamlit UI
    

subheader = st.empty()

# if 'records_keys' not in st.session_state or st.session_state.get('records_keys', {}) == []:
#     st.session_state.records_keys = list(st.session_state.get('records', {}).keys())

# records_keys = st.session_state.get('records_keys', [])
# if len(records_keys) > 0:
if 'records_map' in st.session_state:
    if 'regions_idx' not in st.session_state:
        st.session_state.regions_idx = 0

    # Sidebar for region navigation
    st.sidebar.markdown("### Select Region to Visualize")
    # activate the variable when pressing inter button
    region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None, key="region", help="Enter the region in the format: chr:start-end")
    #_level = st.sidebar.slider('Zoom Level', min_value=1, max_value=100, value=100)
    display_option = st.sidebar.radio("Select Display Type", 
                                ("Sequence with Highlighted Motifs", "Bars"))

    col1, _, col2 = st.columns([1,2, 1])  # Adjust the ratio [1, 1] to control spacing between buttons

    # Place the "Previous region" and "Next region" buttons in these columns
    with col1:
        if st.button("Previous region"):
            region = None
            st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

    with col2:
        if st.button("Next region"):
            region = None
            st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(st.session_state.records_map) - 1)
        
    if region:

        try:
            chr_input, start_end_input = region.split(':')
            start_input, end_input = map(int, start_end_input.split('-'))
            start_input-=1
            end_input-=1

            input_region = f"{chr_input}:{start_input}-{end_input}"
            record_key = st.session_state.records[input_region]

        except:
            st.sidebar.info("Invalid region format, showing the first record")
            record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]
    else:
        st.write(st.session_state.records_map[st.session_state.regions_idx])
        record_key = st.session_state.records[st.session_state.records_map[st.session_state.regions_idx]]


    record = parse_record(vcf_file_path, record_key)
    if len(record["motif_ids_h1"]) == 0 and len(record["motif_ids_h2"]) == 0:
        st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
        st.stop()
        
    subheader = st.subheader(f"TR-region: {record_key}")
    left_column, right_column = st.columns([4, 1])
    # define the motif colors

    motif_colors = get_color_palette(len(record['motifs']))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}

    if display_option == "Sequence with Highlighted Motifs":
        visulize_TR_with_dynamic_sequence(record, left_column, right_column,motif_colors)
    else:
        display_motifs_with_bars(record, left_column, right_column,motif_colors)

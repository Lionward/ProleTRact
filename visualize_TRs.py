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


# Function to parse the VCF record
def parse_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    
    records = {}
    

    for  rec in vcf.fetch():
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
    # iterate over the motifs and set them to 0 if they are not in the motif_count
    for motif in motif_names:
        if motif not in found_motifs_h1:
            motif_count_h1[motif_names.index(motif)] = 0
        if motif not in found_motifs_h2:
            motif_count_h2[motif_names.index(motif)] = 0
    
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
        st.write(label)
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
    motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][1])
   
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
            f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
            f' Motif {idx}: {motif}</div>', unsafe_allow_html=True)
    # show the inturruptions as well as the gray catagory
    st.markdown(
        f'<div id="legend-motif-interruption" class="legend-item motif-interruption" style="background-color:#FF0000;color:black;padding:5px;margin-bottom:10px;border-radius:5px;'
        f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
        f' Interruption</div>', unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)


# Streamlit UI

st.sidebar.title("Tandem Repeat Visualization")


vcf_file = st.sidebar.file_uploader("Upload VCF file", type=["vcf", "gz"])

if 'records' not in st.session_state and vcf_file:
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        tmp_file.write(vcf_file.read())
        tmp_file_path = tmp_file.name
    st.session_state.records = parse_vcf(tmp_file_path)

records = st.session_state.get('records', {})
# definde a container for the subheader

# Example usage in Streamlit UI
subheader = st.empty()




if records:
    # clear all the previous visualizations by clearing the page

    st.sidebar.markdown(f"Total number of Records: {len(records)}")
    records_keys = list(records.keys())

    if 'regions_idx' not in st.session_state:
        st.session_state.regions_idx = 0

    # Sidebar for region navigation
    st.sidebar.markdown("### Select Region to Visualize")
    region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None)
    #_level = st.sidebar.slider('Zoom Level', min_value=1, max_value=100, value=100)
    display_option = st.sidebar.radio("Select Display Type", 
                                  ("Sequence with Highlighted Motifs", "Bars"))

    col1, col2 = st.sidebar.columns([1, 1])  # Adjust the ratio [1, 1] to control spacing between buttons

    # Place the "Previous region" and "Next region" buttons in these columns
    with col1:
        if st.button("Previous region"):
            region = None
            st.session_state.regions_idx = max(st.session_state.regions_idx - 1, 0)

    with col2:
        if st.button("Next region"):
            region = None
            st.session_state.regions_idx = min(st.session_state.regions_idx + 1, len(records_keys) - 1)
        
    if region:
        try:
            chr_input, start_end_input = region.split(':')
            start_input, end_input = map(int, start_end_input.split('-'))
            record_key = f"{chr_input}:{start_input}-{end_input}"
        except:
            st.sidebar.info("Invalid region format, showing the first record")
            record_key = records_keys[st.session_state.regions_idx]
    else:
        record_key = records_keys[st.session_state.regions_idx]


    record = records[record_key]
    subheader = st.subheader(f"TR-region: {record_key}")
    left_column, right_column = st.columns([4, 1])
    # define the motif colors

    motif_colors = get_color_palette(len(record['motifs']))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}

    if display_option == "Sequence with Highlighted Motifs":
        visulize_TR_with_dynamic_sequence(record, left_column, right_column,motif_colors)
    else:
        display_motifs_with_bars(record, left_column, right_column,motif_colors)
    #
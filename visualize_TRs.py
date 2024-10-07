import streamlit as st
import pysam
import re
import matplotlib.pyplot as plt
import pandas as pd
import altair as alt
st.set_page_config(layout="wide")
import os
import pysam
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
import plotly.graph_objects as go
import plotly.express as px
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
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
    if isinstance(vcf_file, str):
        vcf = pysam.VariantFile(vcf_file)
    else:
        vcf = vcf_file
    rec = vcf.fetch(region=region)
    # get the record with the id
    for rec in vcf.fetch(region=region):
        break
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
            alt_allele2 = ''
    if ids_h1 is None:
        ids_h1 = []
    if ids_h2 is None:
        ids_h2 = []
    CN_H1 = rec.info['CN_H1']
    CN_H2 = rec.info['CN_H2']
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
            'motif_ids_h1': ids_h1,
            'motif_ids_h2': ids_h2,
            'motif_ids_ref': rec.info['MOTIF_IDs_REF'],
            'ref_CN': ref_CN,
            'CN_H1': CN_H1,
            'CN_H2': CN_H2,
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

    if record['alt_allele2'] != '':
        CN2_col.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 2 Total copy number:</strong> {str(record['spans'][2]).count('-')}
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
            display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names)
            if record['alt_allele2'] != '':
                display_motifs_as_bars("Allel2",motif_colors, record['motif_ids_h2'], record['spans'][2], record['alt_allele2'], motif_names)
        with tab1:
            display_motifs_as_bars("Allel1",motif_colors, record['motif_ids_h1'], record['spans'][1], record['alt_allele1'], motif_names)
            plot_container_h1 = st.empty()
            with plot_container_h1:
                if show_comparison == False:
                    plot_motif_bar(motif_count_h1, motif_names, motif_colors)
            
            if record['alt_allele2'] != '':
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
    motif_count_ref = count_motifs(record['motif_ids_ref'], record['spans'][0])
    found_motifs_ref = list(motif_count_ref.keys())
    found_motifs_ref = [motif_names[int(m)] for m in found_motifs_ref]
    motif_count_h1 = count_motifs(record['motif_ids_h1'], record['spans'][1])
    found_motifs_h1 = list(motif_count_h1.keys())
    found_motifs_h1 = [motif_names[int(m)] for m in found_motifs_h1]
    motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][2])
    found_motifs_h2 = list(motif_count_h2.keys())
    found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
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
  
    for idx, sequence in enumerate(sequences):
        sequence_name = sequence['name']
        motif_ids_seq = motif_ids[idx]
        spans = spans_list[idx]
        ranges = parse_motif_range(spans)
        sequence_length = len(sequence['sequence'])
        previous_end = 0
 
        for i, (start, end) in enumerate(ranges):
            motif = motif_ids_seq[i]
            color = motif_colors[int(motif)]

            # Handle interruption between motifs
            if start > previous_end:
                data.append({
                    'sample': sequence_name,
                    'Start': previous_end,
                    'End': start,
                    'Motif': 'Interruption',
                    'Color': '#FF0000',
                    'Sequence': sequence['sequence'][previous_end:start],
                })

            # Add motif data
            data.append({
                'sample': sequence_name,
                'Start': start,
                'End': end + 1,  # Altair works with non-inclusive end
                'Motif': motif_names[int(motif)],
                'Color': color,
                'Sequence': sequence['sequence'][start:end+1],
            })

            previous_end = end + 1

        # Add interruption after the last motif if any
        if previous_end < sequence_length:
            data.append({
                'sample': sequence_name,
                'Start': previous_end,
                'End': sequence_length,
                'Motif': 'Interruption',
                'Color': '#FF0000',
                'Sequence': sequence['sequence'][previous_end:],
            })

    return pd.DataFrame(data)


def plot_Cohort_results(cohort_records):
    sequences = []
    span_list = []
    motif_ids_list = []
    for key in cohort_records.keys():
        sequences.append({'name': f'{key}_alle1', 'sequence': cohort_records[key]['alt_allele1']})
        span_list.append(cohort_records[key]['spans'][1])
        motif_ids_list.append(cohort_records[key]['motif_ids_h1'])
        if cohort_records[key]['alt_allele2'] != '':
            sequences.append({'name': f'{key}_alle2', 'sequence': cohort_records[key]['alt_allele2']})
            span_list.append(cohort_records[key]['spans'][2])
            motif_ids_list.append(cohort_records[key]['motif_ids_h2'])

    motif_names = cohort_records[list(cohort_records.keys())[0]]['motifs']
    record = cohort_records[list(cohort_records.keys())[0]]
    motif_colors, df = stack_plot(record, motif_names, sequences, span_list, motif_ids_list)

    # Filterung der Daten
    figure = go.Figure()

    # Einzigartige Proben erhalten
    unique_samples = df['sample'].unique()
    # Unterbrechungen entfernen
    unique_samples = [sample for sample in unique_samples if sample != "Interruption"]

def plot_HGSVC_VS_allele(record, hgsvc_records, motif_names):
    sequences = []
    span_list = []
    motif_ids_list = []

    sequences.append({'name': "Ref", 'sequence': record['ref_allele']})
    span_list.append(record['spans'][0])
    motif_ids_list.append(record['motif_ids_ref'])
    sequences.append({'name': "Allel1", 'sequence':  record['alt_allele1']})
    span_list.append(record['spans'][1])
    motif_ids_list.append(record['motif_ids_h1'])
    if record['alt_allele2'] != '':
        sequences.append({'name': "Allel2", 'sequence': record['alt_allele2']})
        span_list.append(record['spans'][2])
        motif_ids_list.append(record['motif_ids_h2'])

    for key in hgsvc_records.keys():
        if hgsvc_records[key]['alt_allele'] == '':
            continue
        sequences.append({'name': key, 'sequence': hgsvc_records[key]['alt_allele']})
        span_list.append(hgsvc_records[key]['spans'][1])
        motif_ids_list.append(hgsvc_records[key]['motif_ids_h'])

    motif_colors, df = stack_plot(record, motif_names, sequences, span_list, motif_ids_list)


    # Filterung der Daten
    figure = go.Figure()

    # Einzigartige Proben erhalten
    unique_samples = df['sample'].unique()
    # Unterbrechungen entfernen
    unique_samples = [sample for sample in unique_samples if sample != "Interruption"]

    # Scatter-Plots für jede Probe und jedes Motiv hinzufügen
    ref_data = []
    allele_data = []

    # Scatter-Plots für jede Probe und jedes Motiv hinzufügen
    for sample in unique_samples:
        sample_df = df[df['sample'] == sample]
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
    figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color='gray', size=20), name='HGSVC'))
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
    figure.update_yaxes(range=[0, df['sample'].value_counts().max()])
    
    # Plotly-Figur in Streamlit anzeigen
    st.plotly_chart(figure, use_container_width=True)
    #Motif Overlap Radial Chart
    # Pivot tables for HGSVC and sample data
    pivot_hgsvc = pd.pivot_table(df[df['sample'] == 'HGSVC'], index='Motif', columns='sample', values='Length', aggfunc='count', fill_value=0)
    pivot_sample = pd.pivot_table(df[df['sample'] != 'HGSVC'], index='Motif', columns='sample', values='Length', aggfunc='count', fill_value=0)

    # Convert pivot tables to long format for Altair
    pivot_hgsvc_long = pivot_hgsvc.reset_index().melt(id_vars='Motif', var_name='Sample', value_name='Count')
    pivot_sample_long = pivot_sample.reset_index().melt(id_vars='Motif', var_name='Sample', value_name='Count')

    # Combine the data for comparison
    combined_data = pd.concat([pivot_hgsvc_long, pivot_sample_long])

    # Create Altair heatmap
    heatmap = alt.Chart(combined_data).mark_rect().encode(
        x=alt.X('Sample:N', title='Sample'),
        y=alt.Y('Motif:N', title='Motif'),
        color=alt.Color('Count:Q', scale=alt.Scale(scheme='reds'), title='Count'),
        tooltip=['Sample', 'Motif', 'Count']
    ).properties(
        width=600,
        height=400,
        title='Motif Occurrences Heatmap'
    )

    # Display heatmap in Streamlit
    st.altair_chart(heatmap, use_container_width=True)

def stack_plot(record, motif_names, sequences, span_list, motif_ids_list):
    motif_colors = get_color_palette(len(record['motifs']))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
    # plot the min and max copy number in a nice way 

  
           
    df = create_motif_dataframe(sequences, motif_colors, motif_ids_list, span_list, motif_names)
    
    df['Length'] = df['End'] - df['Start'] 
            
            # sort the data frame by start and end per sequence

    df['Order'] = df.index  # Use the index to maintain the original order from the DataFrame
            # Create the Altair chart with explicit order encoding
    default_hight = 600 
    chart_height = max(default_hight, len(sequences) * 7)
    
    df['sample'] = df['sample'].apply(lambda x: x.replace("_pathogenic", ""))
    # sort the samples by name 
    df['sample'] = pd.Categorical(df['sample'], categories=sorted(df['sample'].unique()), ordered=True)
    min_copy_number = df.groupby('sample')['Length'].sum().min()
    max_copy_number = df.groupby('sample')['Length'].sum().max()

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
    chart = alt.Chart(df).mark_bar().encode(
                y=alt.Y(
                    'sample', 
                    sort=None, 
                    axis=alt.Axis(labelOverlap=False, ticks=False)  # Ensure all labels are shown, and avoid overlapping
                ),
                x=alt.X('Length', title='Length', stack='zero'),  # Stack motifs without sorting them
                color=alt.Color('Motif', scale=alt.Scale(domain=list(motif_names) + ['Interruption'], range=list(motif_colors.values()) + ['#FF0000'])),
                order=alt.Order('Order', sort='ascending'),  # Explicitly order the bars by the 'Order' column
                tooltip=['sample', 'Motif', 'Start', 'End', 'Sequence']
            ).properties(
                width=800,
                height=chart_height,  # Use the dynamically calculated chart height
                title="Motif Occurrences"
            ).configure_axis(
                labelFontSize=10,  # Adjust label font size if necessary to fit more labels
                titleFontSize=12
            )
    if chart_height > default_hight:
                # remove y axis labels
        chart = chart.configure_axisY(labelFontSize=0)
    
            # Display the chart in Streamlit
    st.altair_chart(chart, use_container_width=True)
    return motif_colors,df

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
    motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][2])
    found_motifs_h2 = list(motif_count_h2.keys())
    found_motifs_h2 = [motif_names[int(m)] for m in found_motifs_h2]
    motif_count_h1 = {int(k): v for k, v in motif_count_h1.items()}
    motif_count_h2 = {int(k): v for k, v in motif_count_h2.items()}
    total_copy_number_h1 = str(record['spans'][1]).count('-')
    CN1_col.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 1 Total copy number:</strong> {total_copy_number_h1}
            </div>
        """, unsafe_allow_html=True)
    if record['alt_allele2'] != '':
            total_copy_number_h2 = str(record['spans'][2]).count('-')
            CN2_col.markdown(f"""
                <div style="font-size: 20px; color: #FF5733;">
                    <strong>Allele 2 Total copy number:</strong> {total_copy_number_h2}
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
            display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names)
            if record['alt_allele2'] != '':
                alt_allele2 = record['alt_allele2'] #if record['alt_allele2'] != "." else record['ref_allele']
                display_dynamic_sequence_with_highlighted_motifs("Allel2",alt_allele2, record['motif_ids_h2'], record['spans'][2], motif_colors, motif_names)
        with tab1:
        # Render the scrollable sequence with highlighted motifs for allele 1
            alt_allele1 = record['alt_allele1']
            display_dynamic_sequence_with_highlighted_motifs("Allel1",alt_allele1, record['motif_ids_h1'], record['spans'][1], motif_colors, motif_names)

            # Create an empty container for the plot to refresh
            plot_container_h1 = st.empty()
            with plot_container_h1:
                if show_comparison == False:
                    plot_motif_bar(motif_count_h1, motif_names, motif_colors)

            if record['alt_allele2'] != '':

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
    <div class='sidebar-container'>
        <h1>Tandem Repeat Visualization</h1>
        <hr>
        <p>Visualize and analyze tandem repeats in your genomic data.</p>
    </div>
""", unsafe_allow_html=True)

# Container for file upload make it smaller
st.sidebar.markdown("""
    <div class='upload-container'>
        <h3>Upload VCF File</h3>
        <p>Please enter the path of your VCF file to get started.</p>
        
    </div>
""", unsafe_allow_html=True)


path_changed = False
old_vcf_file_path = st.session_state.get('vcf_file_path', None)


st.session_state.analysis_mode = st.sidebar.radio("Select the type of analysis", ("indivisual sample", "Cohort"))

if st.session_state.analysis_mode == "indivisual sample":
    st.sidebar.text_input("", value=None, key="vcf_file_path")
    vcf_file_path = st.session_state.get('vcf_file_path', None)
    if vcf_file_path is None:
        st.stop()
    if vcf_file_path != old_vcf_file_path:
        path_changed = True
        st.session_state.vcf_file_path = vcf_file_path
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
                st.session_state.records,st.session_state.records_map = parse_vcf(vcf_file_path)
                st.session_state.hgsvc_path = "/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output/hgsvc/TandemTwist/asm/"
                #samples_keys = list(st.session_state.get('hgsvc_pop_records', {}).keys())
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
                                    ("Sequence with Highlighted Motifs", "Bars"))

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
            
        if region:

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


        record = parse_record(vcf_file_path, record_key)
        hgsvc_records = get_results_hgsvc_pop(record_key, st.session_state.files ,st.session_state.file_paths)
        #st.session_state.hgsvc_pop_records, = parse_hgsvc_pop(vcf_status,region)
        if len(record["motif_ids_h1"]) == 0 and len(record["motif_ids_h2"]) == 0:
            st.warning(f"No motifs found in the region: {st.session_state.records_map[st.session_state.regions_idx]}")
            st.stop()
            
        middel.markdown(f"""
            <div style="font-size: 18px; color: #90EE90; margin-bottom: 5px; text-align: center; 
                        padding: 5px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1); 
                        background-color: #333; display: inline-block;">
                <strong>Tandem Repeat Region: {record['chr']}:{record['pos']-1}-{record['stop']-1}</strong>
            </div>
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
    

else:
    st.sidebar.text_input("Enter the path to the cohort results", value=None, key="cohort_path")
    st.session_state.path_to_cohort = st.session_state.get('cohort_path', None)
    if st.session_state.path_to_cohort is None:
        st.stop()
    st.session_state.cohort_file_paths = [f for f in os.listdir(st.session_state.path_to_cohort) if f.endswith('CCS.vcf.gz')]
    st.session_state.cohort_files = [load_vcf(st.session_state.path_to_cohort + f) for f in st.session_state.cohort_file_paths]
    st.session_state.cohorts_records_map = get_records_info(st.session_state.path_to_cohort+ st.session_state.cohort_file_paths[0])
    col1, middel, col2 = st.columns([1.5,3, 1])  # Adjust the ratio [1, 1] to control spacing between buttons
    # Place the "Previous region" and "Next region" buttons in these columns
    if 'cohorts_records_map' in st.session_state:
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

        region = st.session_state.cohorts_records_map[st.session_state.regions_idx]
        st.session_state.cohort_results = get_results_cohort(region, st.session_state.cohort_files, st.session_state.cohort_file_paths)
        if 'cohort_results' in st.session_state:   
            region = st.session_state.regions_idx
            plot_Cohort_results(st.session_state.cohort_results)
        else:
            st.stop()


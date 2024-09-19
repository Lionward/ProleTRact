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

# import the browser

# Function to parse the VCF record
def parse_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    
    records = {}
    
    progress = st.progress(0)
    i = 0
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
     
        record = {
            'chr': rec.chrom,
            'pos': rec.pos,
            'motifs': rec.info['MOTIFS'],
            'motif_ids_h1': ids_h1,
            'motif_ids_h2': ids_h2,
            'spans': rec.samples[0]['SP'],
            'ref_allele': ref_allele,
            'alt_allele1': alt_allele1,
            'alt_allele2': alt_allele2
        }
        records[rec.id] = record
        progress.progress((i+1)/1243954)
        i += 1
    progress.empty()

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

# Function to display motifs relative to the sequence, with spans aligned horizontally
def display_motifs_with_bars(motif_ids, spans, motif_colors, sequence_length, motif_names):
    ranges = parse_motif_range(spans)
    if not isinstance(motif_names, list):
        motif_names = [motif_names]
    # make the border thicker
    bar_container = "<div style='width:100%; position: absolute; height: 20px; border:2px solid black; border-radius: 8px;'>"
    previous_end = 0
    gap = 0.05

    # Loop through each motif in the sequence
    for idx, (start, end) in enumerate(ranges):
        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        span_length = end - start + 1
        
        if start >= 0 and end <= sequence_length:
            if start > previous_end:
                interruption_width = (start - previous_end) / sequence_length * 100
                interruption_start = previous_end / sequence_length * 100
                bar_container += (
                    f"<div style='position:absolute; background-color:#D3D3D3; left:{interruption_start}%; "
                    f"width:{interruption_width}%; height:18px;top:-1px; border-radius:6px  ;border:1px solid black;'></div>"
                )
            
            relative_width = (span_length / sequence_length) * 100 - gap
            relative_start = (start / sequence_length) * 100
            
            # Assign each motif a unique class for highlight
            bar_container += (
                f"<div id='legend-motif-{int(motif)}' class='motif-sequence motif-{idx}' style='position:absolute; background-color:{color}; left:{relative_start}%; "
                f"width:{relative_width}%; height:18px; top:-1px; border-radius:6px  ;border:1px solid black; "
                f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2);"
                f"cursor: pointer;' title='Motif: {motif_names[int(motif)]}'>"
                f"</div>"
            )
            previous_end = end + 1
    
    if previous_end < sequence_length:
        interruption_width = (sequence_length - previous_end) / sequence_length * 100
        interruption_start = previous_end / sequence_length * 100
        bar_container += (
            f"<div style='position:absolute; background-color:#D3D3D3; left:{interruption_start}%; "
            f"width:{interruption_width}%; height:18px; top:21px; border-radius:8px;'></div>"
        )
    
    bar_container += "</div>"
    
    st.markdown(bar_container, unsafe_allow_html=True)

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
    
    # Loop through each motif in the list and count occurrences
    for idx, motif in enumerate(motif_ids):
        if motif in motif_count:
            motif_count[motif] += 1
        else:
            motif_count[motif] = 1
    
    return motif_count

# Function to plot motif occurrences as a bar chart
def plot_motif_bar(motif_count, motif_names, motif_colors=None):
    # Convert motif indices to names
    motif_labels = [motif_names[int(m)] for m in motif_count.keys()]
    
    # Create a DataFrame for plotting
    data = {
        'Motif': motif_labels,
        'Count': list(motif_count.values())
    }
    df = pd.DataFrame(data)
    colors = [motif_colors[int(motif)] for motif in motif_count.keys()]
    
    # Plot the motif occurrences using Plotly
    fig = px.bar(
        df,
        x='Motif',
        y='Count',
        color='Motif',
        color_discrete_sequence=colors,
        title="Motif Occurrences",
        labels={'Count': 'Motif Count'},
        width=800,
        height=400
    )
    
    # Customize the layout
    fig.update_layout(
        xaxis_title="Motif",
        yaxis_title="Count",
        xaxis_tickangle=-45,
        bargap=0.2,
        height=200,
        margin=dict(t=50, b=50)
    )
    
    # Annotate the bars with the count values
    fig.update_traces(texttemplate='%{y}', textposition='inside')
    # remove the legend
    fig.update_layout(showlegend=False)
    # remove the x-axis ticks
    fig.update_xaxes(showticklabels=False)
    # Display the plot in Streamlit
    st.plotly_chart(fig)

# def plot_motif_bar(motif_count, motif_names, motif_colors=None):
#     # Convert motif indices to names
#     motif_labels = [motif_names[int(m)] for m in motif_count.keys()]


    
#     # Create a DataFrame for plotting
#     data = {
#         'Motif': motif_labels,
#         'Count': list(motif_count.values())
#     }
#     df = pd.DataFrame(data)
#     colors = []
#     for idx, motif in enumerate(motif_count.keys()):
#         colors.append(motif_colors[int(motif)])
    
#     # Plot the motif occurrences
#     fig, ax = plt.subplots(figsize=(8, 2))  # Make the plot smaller
#     sns.barplot(x='Motif', y='Count', data=df, palette=colors, ax=ax)  # Swap x and y
#     ax.set_title("Motif Occurrences")
#     ax.set_xlabel("Motif")
#     ax.set_ylabel("Count")
#     # addjust the width of the bars and locate them in the middle

#     for bar in ax.patches:
#         bar.set_width(0.5)
#         bar.set_x(bar.get_x() + 0.15)

#     plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
    
#     # Annotate the bars with the count values
#     for p in ax.patches:
#         ax.annotate(format(p.get_height(), '.1f'), 
#                     (p.get_x() + p.get_width() / 2., p.get_height()-max(df['Count'])/4), 
#                     ha = 'center', va = 'center', 
#                     xytext = (0, 9), 
#                     textcoords = 'offset points')
    
#     plt.tight_layout()
#     # move the plot a bit down
#     plt.subplots_adjust(top=1.1)

#     st.pyplot(fig)


# Adjust the visualization to zoom into a certain portion of the sequence
def display_zoomed_motifs(motif_ids, spans, motif_colors, sequence, motif_names, zoom_level):
    ranges = parse_motif_range(spans)
    sequence_length = len(sequence)
    zoom_length = int(sequence_length * (zoom_level / 100))  # Adjust visible length by zoom level
    start_pos = 0  # This can be adjusted dynamically later for scrolling or dragging


    if  isinstance(motif_names, tuple):
        motif_names = list(motif_names)
    elif not isinstance(motif_names, list):
        motif_names = [motif_names]
        
    # Only render motifs within the zoomed-in window
    bar_container = "<div style='width:100%; position: relative; height: 20px; border:2px solid black; border-radius: 8px;'>"
    previous_end = 0
    gap = 0.05
    #2st.write(f"sequence : {sequence}")
    for idx, (start, end) in enumerate(ranges):
        #st.write(f"start: {start} end: {end}")
    
        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        span_length = end - start +1
        
        if start >= start_pos and end <= start_pos + zoom_length:
            if start > previous_end:
                interruption_width = (start - previous_end) / zoom_length * 100
                interruption_start = previous_end / zoom_length * 100
                bar_container += (
                    f"<div style='position:absolute; background-color:#FF0000; left:{interruption_start}%; "
                    f"width:{interruption_width}%; height:18px;top:-1px; border-radius:6px;"
                    f"cursor: pointer;' title='Inturrutoptions: {sequence[previous_end:start]}'>"
                    f"</div>"
                )
                
            
            relative_width = (span_length / zoom_length) * 100 - gap
            relative_start = (start - start_pos) / zoom_length * 100
            
            bar_container += (
                f"<div id='legend-motif-{int(motif)}' class='motif-sequence motif-{idx}' style='position:absolute; background-color:{color}; left:{relative_start}%; "
                f"width:{relative_width}%; height:18px; top:-1px; border-radius:6px  ;border:1px solid black; "
                f"transition: all 0.3s ease; box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2);"
                f"cursor: pointer;' title='Motif: {motif_names[int(motif)]}'>"
                f"</div>"
            )
            previous_end = end +1

    if previous_end < zoom_length:
        interruption_width = (zoom_length - previous_end) / zoom_length * 100
        interruption_start = previous_end / zoom_length * 100
        bar_container += (
            f"<div style='position:absolute; background-color:#D3D3D3; left:{interruption_start}%; "
            f"width:{interruption_width}%; height:18px;top:-1px; border-radius:6px;'></div>"
        )
    
    bar_container += "</div>"
    st.markdown(bar_container, unsafe_allow_html=True)

# Function to display the motif legend on the right side
def display_motif_legend(motifs, motif_colors, right_column):
    with right_column:
        st.markdown("### Motif Legend")
        st.markdown('<div style="max-height:400px; overflow-y:scroll;">', unsafe_allow_html=True)  # Scrollable container
        if  isinstance(motifs, tuple):
            motifs = list(motifs)
        elif not isinstance(motifs, list):
            motifs = [motifs]

        for idx, motif in enumerate(motifs):
            color = motif_colors[idx]
            # Assign each legend item a unique class for hover effect
            right_column.markdown(
                f'<div id="legend-motif-{idx}" class="legend-item motif-{idx}" style="background-color:{color};color:white;padding:5px;margin-bottom:10px;border-radius:5px;'
                f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
                f' Motif {idx}: {motif}</div>', unsafe_allow_html=True)
        # show the inturruptions as well as the gray catagory
        right_column.markdown(
            f'<div id="legend-motif-interruption" class="legend-item motif-interruption" style="background-color:#FF0000;color:black;padding:5px;margin-bottom:10px;border-radius:5px;'
            f' text-align:center;font-weight:bold;font-size:12px;border:2px solid #000;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.3);">'
            f' Interruption</div>', unsafe_allow_html=True)
        
        right_column.markdown('</div>', unsafe_allow_html=True)


st.markdown("""
    <style>
    /* CSS to highlight corresponding motifs on hover */
    .highlight { 
        transform: scale(1.2); 
        box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.5);
        z-index: 10;
    }

    /* CSS to make non-hovered motifs more transparent */
    .transparent {
        opacity: 0.3;
        transition: opacity 0.3s ease;
    }
    </style>

    <script>
    // Function to handle hover effects
    function attachHoverEffects() {
        const motifs = document.querySelectorAll('.legend-item');

        motifs.forEach(motif => {
            const motifClass = motif.classList[1];  // Get the class that links the motif

            // Add mouseover event to highlight corresponding sequence motifs and dim others
            motif.addEventListener('mouseover', () => {
                const allMotifs = document.querySelectorAll('.motif-sequence');
                const sequenceMotifs = document.querySelectorAll(`.${motifClass}`);

                // Make all motifs more transparent
                allMotifs.forEach(m => m.classList.add('transparent'));

                // Highlight only the corresponding motifs
                sequenceMotifs.forEach(m => {
                    m.classList.remove('transparent');
                    m.classList.add('highlight');
                });
            });

            // Add mouseout event to remove the highlight and restore opacity
            motif.addEventListener('mouseout', () => {
                const allMotifs = document.querySelectorAll('.motif-sequence');
                const sequenceMotifs = document.querySelectorAll(`.${motifClass}`);

                // Remove transparency from all motifs
                allMotifs.forEach(m => m.classList.remove('transparent'));

                // Remove highlight from corresponding motifs
                sequenceMotifs.forEach(m => m.classList.remove('highlight'));
            });
        });
    }

    // Use MutationObserver to detect when new elements are added to the DOM
    const observer = new MutationObserver((mutationsList, observer) => {
        for(const mutation of mutationsList) {
            if (mutation.type === 'childList' && mutation.addedNodes.length > 0) {
                attachHoverEffects();  // Attach hover effects once new elements are rendered
            }
        }
    });

    // Observe the entire document for changes
    observer.observe(document.body, { childList: true, subtree: true });
    </script>
""", unsafe_allow_html=True)


# Streamlit UI

st.sidebar.title("Tandem Repeat Visualization")



# Upload VCF file in the sidebar
vcf_file = st.sidebar.file_uploader("Upload VCF file", type=["vcf", "gz"])
# load the index file as well

def visulize_TR(zoom_level, left_column, right_column, record_key, record):

    motif_names = record['motifs']
    if  isinstance(motif_names, tuple):
        motif_names = list(motif_names)
    elif not isinstance(motif_names, list):
        motif_names = [motif_names]

    
    motif_colors = get_color_palette(len(motif_names))
    motif_colors = {idx: color for idx, color in enumerate(motif_colors)}
      # These are used in the labels
    motif_count_h1 = count_motifs(record['motif_ids_h1'], record['spans'][0])
    motif_count_h2 = count_motifs(record['motif_ids_h2'], record['spans'][1])
            # Display motif legend in the right column
    display_motif_legend(motif_names, motif_colors, right_column)

    with left_column:
        st.subheader(f"TR-region: {record_key}")

        # Allele 1 Display (Motifs displayed as bars relative to sequence length)
        total_copy_number_h1 = str(record['spans'][0]).count('-')

        # Verwenden von HTML und CSS, um die Textfarbe zu ändern
        st.markdown(f"""
            <div style="font-size: 20px; color: #FF5733;">
                <strong>Allele 1 Total copy number:</strong> {total_copy_number_h1}
            </div>
        """, unsafe_allow_html=True)
        display_zoomed_motifs(record['motif_ids_h1'], record['spans'][0], motif_colors, record['alt_allele1'], motif_names, zoom_level)

        plot_motif_bar(motif_count_h1, motif_names, motif_colors)
                # Allele 2 Display (if available)
        if record['alt_allele2'] != '':
            total_copy_number_h2 = str(record['spans'][1]).count('-')
            st.markdown(f"""
                <div style="font-size: 20px; color: #FF5733;">
                    <strong>Allele 2 Total copy number:</strong> {total_copy_number_h2}
                </div>
            """, unsafe_allow_html=True)

            display_zoomed_motifs(record['motif_ids_h2'], record['spans'][1], motif_colors, record['alt_allele2'], motif_names, zoom_level)
                    # make distance between the plots
            st.markdown('<vr>', unsafe_allow_html=True)
            plot_motif_bar(motif_count_h2, motif_names, motif_colors)


def get_color_palette(n):
        # Using the 'viridis' colormap from matplotlib
    cmap = plt.get_cmap('tab20')  # You can try 'viridis', 'plasma', 'Set1', etc.
    colors = [cmap(i) for i in range(n)]
    return ['#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in colors]
# Check if the VCF file has been parsed before
if 'records' not in st.session_state:
    if vcf_file:
        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(vcf_file.read())
            tmp_file_path = tmp_file.name
        # Parse the VCF file and store it in session state
        st.session_state.records = parse_vcf(tmp_file_path)


records = st.session_state.get('records', {})

# Further processing and visualization
if records:
    # Your visualization code continues here
    # For example, display summary statistics
    st.sidebar.markdown(f"Total number of Records: {len(records)}")
    records_keys = list(records.keys())
    # check if  st.session_state.regions_idx is bigger than the number of records
    if 'regions_idx' in st.session_state:
        if st.session_state.regions_idx < 0:
            st.session_state.regions_idx = 0
        if st.session_state.regions_idx >= len(records_keys):
            st.session_state.regions_idx = len(records_keys) - 1



    # clear all the previous values
    
    # Sidebar inputs for selective region visualization
    st.sidebar.markdown("### Select Region to Visualize")
    region = st.sidebar.text_input("TR region (e.g., chr1:1000-2000)", value=None)
    zoom_level = st.sidebar.slider('Zoom Level', min_value=1, max_value=100, value=100)
    
    # Initialize regions_idx in session state if not already done
    if 'regions_idx' not in st.session_state:
        st.session_state.regions_idx = 0
    if st.sidebar.button("Next region"):
        st.session_state.regions_idx += 1  
    if st.sidebar.button("Previous region"):
        st.session_state.regions_idx -= 1
    
    # split region into chr, start, end
    if region:
        chr_input, start_end_input = region.split(':')
        start_input, end_input = start_end_input.split('-')
        start_input = int(start_input)
        end_input = int(end_input)
    else:
        chr_input = ''
        start_input = 0
        end_input = 0
    
    # Create two columns: one for visualization, one for legend
    left_column, right_column = st.columns([4, 1])  # Make left column wider than right column

    # Generate color palette for motifs
    if start_input != 0 and end_input != 0:
        record_key = f"{chr_input}:{start_input}-{end_input}"

        if record_key in records:
            st.session_state.regions_idx = records_keys.index(record_key)
            record = records[record_key]

            visulize_TR( zoom_level, left_column, right_column, record_key, record) 
            # Clear the previous plot
            left_column.empty()
            right_column.empty()
        else:
            # Clear the previous plot
            left_column.empty()
            right_column.empty()
            # Show the first record
            st.sidebar.info("Region not found, showing the first record")
            record_key = records_keys[st.session_state.regions_idx]
            record = records[record_key]

            visulize_TR( zoom_level, left_column, right_column, record_key, record)
    else:
        # Show the first record
        if st.session_state.regions_idx < 0 or st.session_state.regions_idx >= len(records_keys):
            st.session_state.regions_idx = 0
            record_key = records_keys[st.session_state.regions_idx]
            record = records[record_key]

            visulize_TR( zoom_level, left_column, right_column, record_key, record)
        else:
            if st.session_state.regions_idx >= len(records_keys):
                st.session_state.regions_idx = len(records_keys) - 1
            record_key = records_keys[st.session_state.regions_idx]
            record = records[record_key]

            visulize_TR( zoom_level, left_column, right_column, record_key, record)



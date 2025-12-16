"""
Plotting functions for visualization.
"""
import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import plotly.graph_objects as go
import plotly.express as px
import re as _re
from sklearn.preprocessing import StandardScaler
from proletract.modules.viz import vis_helper as _vh
from proletract.modules.viz import utils
from proletract.modules.viz import parsers
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd
create_genotype_comparison_matrix = _vh.create_genotype_comparison_matrix


def cluster_plot_samples(df, motif_colors):
        """
        Cluster plot reimplemented using Plotly (static plot).

        Displays clustering of samples by copy number, length, or both.
        Highlights current samples with stars/arrows.
        """

        df = df[(df["Motif"] != "Interruption") & (df["Sample"] != "Interruption")]
        # If no samples available for clustering, show a warning
        if df.empty:
            st.warning("No samples available for clustering.")
            return

        # Get current sample names
        current_sample_names = {"Current Sample", "Allele1", "Allele2", "Allel1", "Allel2", "Reference"}
        # Get current samples
        current_samples = [s for s in df["Sample"].unique() if s in current_sample_names]

        # Compute info per sample (copy number and length)
        def sample_stats(s):
            sub = df[df["Sample"] == s]
            return {
                "Sample": s,
                "Copy Number": len(sub),
                "Length": sub["Length"].sum() if "Length" in sub.columns else 0
            }
        # Create a dataframe of sample stats
        cluster_df = pd.DataFrame([sample_stats(s) for s in df["Sample"].unique()])
        # If less than 2 samples, show a warning
        if len(cluster_df) < 2:
            st.info("Need at least 2 samples for clustering.")
            return

        # Cluster settings (copy number, length, or both)
        cluster_by = st.radio("Cluster by", ["Copy Number", "Length", "Both"], index=2, horizontal=True)

        # Prepare input X (copy number and length)
        if cluster_by == "Both":
            X = cluster_df[["Copy Number", "Length"]].values
            x_col, y_col = "Copy Number", "Length"
            scaler = StandardScaler()
            X_use = scaler.fit_transform(X)
            is_2d = True
            # Add jitter to prevent overlapping points (randomly offset the points)
            rng = np.random.default_rng(seed=42)
            x_range = cluster_df[x_col].max() - cluster_df[x_col].min()
            y_range = cluster_df[y_col].max() - cluster_df[y_col].min()
            jitter_x = rng.normal(0, x_range * 0.01, len(cluster_df))
            # Make y-separation regular so points are clearly vertically stacked above each other 
            sep = y_range * 0.09 if y_range > 0 else 1
            jitter_y = np.arange(len(cluster_df)) * sep
            cluster_df["Copy Number_Jitter"] = cluster_df[x_col] + jitter_x
            cluster_df["Length_Jitter"] = cluster_df[y_col].min() + jitter_y
        else:
            colname = "Copy Number" if cluster_by == "Copy Number" else "Length"
            X = cluster_df[[colname]].values
            x_col = colname
            y_col = None
            X_use = X
            is_2d = False


        min_k = 2
        max_k = min(8, len(cluster_df)) if len(cluster_df) > 2 else 2  # Limit max K for small samples/fast scoring
        inertia_scores = []
        silhouette_scores = []
        candidate_ks = range(min_k, max_k + 1)

        for k in candidate_ks:
            kmeans_test = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = kmeans_test.fit_predict(X_use)
            inertia_scores.append(kmeans_test.inertia_)
            # Silhouette score only valid if n_clusters > 1 and < n_samples
            if k < len(cluster_df):
                try:
                    sil_score = silhouette_score(X_use, labels)
                except Exception:
                    sil_score = float('nan')
            else:
                sil_score = float('nan')
            silhouette_scores.append(sil_score)



        metrics_df = pd.DataFrame({
            'k': list(candidate_ks),
            'Silhouette Score': silhouette_scores,
            'Inertia': inertia_scores
        })


        # Suggest the K with max silhouette, fallback to elbow (lowest inertia after drop)
        if any(pd.notnull(silhouette_scores)):
            try:
                suggested_k = metrics_df.loc[
                    metrics_df['Silhouette Score'] == max([v for v in silhouette_scores if pd.notnull(v)])
                ]['k'].iloc[0]
            except Exception:
                suggested_k = min_k
        else:
            # Use elbow: look for the "elbow" by inertia difference
            diffs = [inertia_scores[i - 1] - inertia_scores[i] for i in range(1, len(inertia_scores))]
            if diffs:
                elbow = diffs.index(max(diffs)) + min_k + 1
                suggested_k = elbow
            else:
                suggested_k = min_k

        num_clusters = st.number_input(
            "Number of Clusters (K)",
            min_value=min_k,
            max_value=max_k,
            value=int(suggested_k),
            step=1,
            help="Choose number of clusters. current cluster number is calculated using silhouette score.",
        )

        kmeans = KMeans(n_clusters=int(num_clusters), random_state=42, n_init=10)
        # Fit KMeans model and predict clusters
        clusters = kmeans.fit_predict(X_use)
        # Shift cluster numbers to start from 1 for display
        cluster_df["Cluster"] = clusters + 1


        # Color palette for clusters
        colors = [
            "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#6A994E", "#7209B7", "#4361EE", "#4CC9F0"
        ]
        # Get unique clusters
        unique_clusters = sorted(cluster_df["Cluster"].unique())
        # Map clusters to colors
        color_map = {c: (colors[i % len(colors)] if c != -1 else "#aaa") for i, c in enumerate(unique_clusters)}
        # Convert cluster numbers to strings
        cluster_df["Cluster_str"] = cluster_df["Cluster"].astype(str)
        # Check if sample is current sample
        cluster_df["Is_Current"] = cluster_df["Sample"].isin(current_samples)

        FONT_SIZE = 24
        FONT_COLOR = '#4B5563'
        MARKER_SIZE_NORMAL = 21
        MARKER_SIZE_CURRENT = 27
        ARROW_COLOR = "#DC2626"
        REF_ARROW_COLOR = "#059669"
        
        # add a red triangle marker for current samples
        def add_arrow(fig, row, x_col, y_col_val, is_2d, color):
            # Compute where to place the triangle below the point for visibility
            TRIANGLE_OFFSET = 0.1  # vertical offset below the point
            TRIANGLE_SIZE = 20      # marker size
            if is_2d:
                point_x = row[f"{x_col}_Jitter"]
                point_y = row[f"{y_col}_Jitter"] + TRIANGLE_OFFSET
            else:
                point_x = row[x_col]
                point_y = row[y_col_val] + TRIANGLE_OFFSET

            # Add a 'triangle-down' marker at the desired location
            fig.add_trace(
                go.Scatter(
                    x=[point_x],
                    y=[point_y],
                    mode="markers",
                    marker=dict(
                        symbol='triangle-down',
                        size=TRIANGLE_SIZE,
                        color=color,
                        line=dict(width=1.5, color="#14532d" if color == REF_ARROW_COLOR else "#aa2020"),
                        opacity=0.99,
                    ),
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

        # Helper function to add cluster traces
        def add_cluster_traces(fig, cluster_df, unique_clusters, color_map, x_col, y_col, is_2d):
            for c in unique_clusters:
                mask = cluster_df["Cluster"] == c
                
                # Separate normal and current samples
                normal_mask = mask & ~cluster_df["Is_Current"]
                current_mask = mask & cluster_df["Is_Current"]
                
                # Add normal cluster points (non-current)
                if normal_mask.any():
                    # Use jittered columns for 2D mode
                    if is_2d:
                        x_data = cluster_df.loc[normal_mask, f"{x_col}_Jitter"]
                        y_data = cluster_df.loc[normal_mask, f"{y_col}_Jitter"]
                    else:
                        x_data = cluster_df.loc[normal_mask, x_col]
                        y_data = cluster_df.loc[normal_mask, "Y_Jitter"]
                    
                    fig.add_trace(go.Scatter(
                        x=x_data,
                        y=y_data,
                        mode="markers",
                        marker=dict(
                            size=MARKER_SIZE_NORMAL,
                            color=color_map[c],
                            opacity=0.85 if is_2d else 0.89,
                            line=dict(width=1.0, color='#eee')
                        ),
                        name=f"Cluster {c}" if c != -1 else "Noise",
                        text=cluster_df.loc[normal_mask, "Sample"],
                        hoverinfo="text",
                        showlegend=True
                    ))
                elif current_mask.any():
                    if is_2d:
                        x_data = cluster_df.loc[current_mask, f"{x_col}_Jitter"].iloc[[0]] 
                        y_data = cluster_df.loc[current_mask, f"{y_col}_Jitter"].iloc[[0]]
                    else:
                        x_data = cluster_df.loc[current_mask, x_col].iloc[[0]]
                        y_data = cluster_df.loc[current_mask, "Y_Jitter"].iloc[[0]]
                    
                    fig.add_trace(go.Scatter(
                        x=x_data,
                        y=y_data,
                        mode="markers",
                        marker=dict(
                            size=MARKER_SIZE_NORMAL,
                            color=color_map[c],
                            opacity=0.01,  # Nearly invisible but ensures legend entry
                            line=dict(width=1.0, color='#eee')
                        ),
                        name=f"Cluster {c}" if c != -1 else "Noise",
                        text=cluster_df.loc[current_mask, "Sample"].iloc[[0]],
                        hoverinfo="skip",  # Skip hover since it's just for legend
                        showlegend=True
                    ))
                
                # Add current samples colored by their cluster
                if current_mask.any():
                    
                    for _, row in cluster_df.loc[current_mask].iterrows():
                        cx = row[f"{x_col}_Jitter"] if is_2d else row[x_col]
                        cy = row[f"{y_col}_Jitter"] if is_2d else row["Y_Jitter"]
                        fig.add_trace(
                            go.Scatter(
                                x=[cx],
                                y=[cy],
                                mode="markers",
                                marker=dict(
                                    size=MARKER_SIZE_CURRENT,
                                    color=color_map[c],
                                    opacity=1.0,
                                    line=dict(width=2.0, color='#0f172a')
                                ),
                                text=[row["Sample"]],
                                hoverinfo="text",
                                showlegend=False
                            )
                        )
                    
                    # Draw arrows for current samples 
                    curr = cluster_df.loc[current_mask].copy()
                    for _, row in curr.iterrows():
                        arrow_color = REF_ARROW_COLOR if row["Sample"] == "Reference" else ARROW_COLOR
                        add_arrow(fig, row, x_col, (y_col if is_2d else "Y_Jitter"), is_2d, arrow_color)
        
        # Initialize figure
        fig = go.Figure()
        
        # Add jitter for 1D case
        if not is_2d:
            rng = np.random.default_rng(seed=42)
            cluster_df["Y_Jitter"] = rng.normal(0, 0.12, len(cluster_df))
        
        # Add all traces
        add_cluster_traces(fig, cluster_df, unique_clusters, color_map, x_col, y_col, is_2d)
        
        # Common axis configuration
        axis_config = dict(
            tickfont=dict(size=FONT_SIZE, color=FONT_COLOR, weight='bold')
        )
        title_font_config = dict(size=FONT_SIZE, color=FONT_COLOR, weight='bold')
        
        # Update layout
        layout_updates = {
            "width": 600,
            "height": 380 if is_2d else 320,
            "xaxis_title": x_col,
            "yaxis_title": y_col if is_2d else "",
            "title": dict(
                text=f"Sample Clustering by {cluster_by} (KMeans)",
                font=dict(size=30, color="#1f2937", family="Arial Black, sans-serif"),
                x=0.5,
                xanchor="center"
            ),
            "legend": dict(
                title="",
                title_font=dict(size=22, color='#374151', family="Arial Black, sans-serif"),
                font=dict(size=18, color='#374151', family="Arial Black, sans-serif"),
                orientation="h",
                y=-1.2,
                yanchor="top",
                x=0.5,
                xanchor="center",
                itemwidth=70
            ),
            "plot_bgcolor": "white",
            "xaxis": {**axis_config, "title": dict(text=x_col, font=title_font_config)},
            "margin": dict(l=60, r=40, t=70, b=110)
        }
        
        if is_2d:
            layout_updates["yaxis"] = {**axis_config, "title": dict(text=y_col, font=title_font_config)}
        else:
            layout_updates["yaxis"] = dict(showline=False, gridcolor='#fafafa', showticklabels=False)
            fig.update_xaxes(showline=True, linewidth=2, linecolor='gray', gridcolor='#ccc', ticks='outside')
        
        fig.update_layout(**layout_updates)

        # Remove modebar to make it non-interactive
        config = {
            "displayModeBar": False,
            "staticPlot": True,
            "displaylogo": True
        }
        st.plotly_chart(fig, use_container_width=True, config=config)

        # Cluster information table
        cluster_info = cluster_df[cluster_df["Cluster"] != -1].groupby("Cluster").agg({
            "Sample": lambda x: list(x),
            "Copy Number": "mean",
            "Length": "mean"
        }).round(2).reset_index()
        
        # Show cluster information table
        if not cluster_info.empty:
            st.html("""
            <style>
                .cluster-section {
                    margin-top: 30px;
                }
                .cluster-title {
                    font-size: 32px;
                    font-weight: 700;
                    color: #1f2937;
                    margin-bottom: 24px;
                    display: flex;
                    align-items: center;
                    gap: 10px;
                }
                .cluster-table {
                    width: 100%;
                    border-collapse: separate;
                    border-spacing: 0;
                    margin-bottom: 20px;
                    background: white;
                    border-radius: 12px;
                    overflow: hidden;
                    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
                }
                .cluster-table thead {
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                }
                .cluster-table th {
                    padding: 20px 24px;
                    text-align: left;
                    font-weight: 700;
                    font-size: 18px;
                    text-transform: uppercase;
                    letter-spacing: 0.5px;
                }
                .cluster-table tbody tr {
                    border-bottom: 1px solid #e5e7eb;
                    transition: background-color 0.2s ease;
                }
                .cluster-table tbody tr:hover {
                    background-color: #f9fafb;
                }
                .cluster-table tbody tr:last-child {
                    border-bottom: none;
                }
                .cluster-table td {
                    padding: 20px 24px;
                    font-size: 17px;
                    color: #374151;
                }
                .cluster-number {
                    display: inline-flex;
                    align-items: center;
                    justify-content: center;
                    width: 48px;
                    height: 48px;
                    border-radius: 10px;
                    font-weight: 700;
                    font-size: 20px;
                    color: white;
                    margin-right: 12px;
                }
                .cluster-number-cell {
                    display: flex;
                    align-items: center;
                }
                .metric-value {
                    font-weight: 700;
                    font-size: 19px;
                    color: #1f2937;
                }
                .metric-label {
                    font-size: 12px;
                    color: #6b7280;
                    margin-right: 8px;
                }
                .samples-cell {
                    display: flex;
                    flex-wrap: wrap;
                    gap: 8px;
                }
                .sample-chip {
                    display: inline-block;
                    background: #f3f4f6;
                    color: #374151;
                    font-size: 15px;
                    font-weight: 600;
                    padding: 7px 14px;
                    border-radius: 6px;
                    border: 1px solid #e5e7eb;
                }
                .noise-box {
                    margin-top: 20px;
                    padding: 18px 24px;
                    background: #fffbeb;
                    border-left: 4px solid #f59e0b;
                    border-radius: 8px;
                    display: flex;
                    align-items: center;
                    gap: 12px;
                }
                .noise-icon {
                    font-size: 24px;
                }
                .noise-text {
                    font-size: 16px;
                    color: #92400e;
                    font-weight: 600;
                }
                .noise-samples-list {
                    color: #78350f;
                    font-weight: 500;
                    font-size: 15px;
                }
            </style>
            """)
            
            st.html('<div class="cluster-section">')
            st.html('<div class="cluster-title">📊 Cluster Summary</div>')
            
            # Create HTML table
            table_html = """
            <table class="cluster-table">
                <thead>
                    <tr>
                        <th>Cluster</th>
                        <th>Copy Number</th>
                        <th>Length</th>
                        <th>Samples</th>
                    </tr>
                </thead>
                <tbody>
            """

            # Show cluster information table
            for _, row in cluster_info.iterrows():
                cluster_num = int(row['Cluster'])
                samples = row['Sample']
                copy_num = row['Copy Number']
                length = row['Length']
                cluster_color = color_map.get(cluster_num, "#667eea")
                
                # Create sample chips
                samples_html = "".join([
                    f'<span class="sample-chip">{sample}</span>' 
                    for sample in samples
                ])
                
                table_html += f"""
                    <tr>
                        <td>
                            <div class="cluster-number-cell">
                                <span class="cluster-number" style="background: {cluster_color};">
                                    {cluster_num}
                                </span>
                            </div>
                        </td>
                        <td>
                            <span class="metric-value">{copy_num:.2f}</span>
                        </td>
                        <td>
                            <span class="metric-value">{length:.2f}</span>
                        </td>
                        <td>
                            <div class="samples-cell">
                                {samples_html}
                            </div>
                        </td>
                    </tr>
                """
            
            table_html += """
                </tbody>
            </table>
            """
            st.html(table_html)
            st.html('</div>')
        
        # Show noise samples
        if (cluster_df["Cluster"] == -1).any():
            noise_samples = cluster_df[cluster_df["Cluster"] == -1]["Sample"].tolist()
            noise_html = f"""
            <div class="noise-box">
                <span class="noise-icon">⚠️</span>
                <div>
                    <div class="noise-text">Outliers (Noise)</div>
                    <div class="noise-samples-list">{", ".join(noise_samples)}</div>
                </div>
            </div>
            """
            st.html(noise_html)


def bar_plot_motif_count(df, region, sort_by="Value"):

        # Prevent overlap
        st.markdown("<div style='height: 40px;'></div>", unsafe_allow_html=True)

        # Extract gene, inheritance, and disease information
        updated_region = region
        pathogenic_threshold = 0
        gene_name = None
        inheritance = None
        disease = None
        
        if updated_region in st.session_state.pathogenic_TRs['region'].unique():
            pathogenic_threshold = st.session_state.pathogenic_TRs.loc[
                (st.session_state.pathogenic_TRs['region'] == updated_region)
            ]['pathogenic_min'].values[0]
            if pathogenic_threshold is None:
                pathogenic_threshold = 0
            gene_name = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['gene'].values[0]
            inheritance = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['inheritance'].values[0]
            disease = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['disease'].values[0]

        df = df[df['Motif'] != "interruption"]
        
        total_copy_number = df.groupby('Sample').size().reset_index(name='Total Copy Number')
        total_copy_number = total_copy_number.sort_values(by='Total Copy Number', ascending=False)
        if sort_by == "Value":
            x_sort = alt.EncodingSortField(field='Total Copy Number', op='sum', order='descending')
        else:
            x_sort = alt.SortField(field='Sample', order='ascending')

        unique_samples = list(total_copy_number['Sample'].apply(lambda x: x.rsplit('_', 1)[0]).unique())
        
        # Updated color palette to match the theme
        color_palette = [
            '#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391',
            '#f6e05e', '#f6ad55', '#fc8181', '#7e9af9', '#c084fc', '#f472b6'
        ]
        color_mapping = {sample: color_palette[i % len(color_palette)] for i, sample in enumerate(unique_samples)}
        color_mapping = {sample: color_mapping[sample.rsplit('_', 1)[0]] for sample in total_copy_number['Sample']}

        # Prepare subtitle text for the chart title
        subtitle_text = ""
        if gene_name or inheritance or disease or updated_region:
            gene_info_parts = []
            if gene_name:
                gene_info_parts.append(f"Gene: {gene_name}")
            if inheritance:
                gene_info_parts.append(f"Inheritance: {inheritance}")
            if disease:
                disease_display = disease
                gene_info_parts.append(f"Disease: {disease_display}")
            if updated_region:
                gene_info_parts.append(f"Region: {updated_region}")
            subtitle_text = " • ".join(gene_info_parts)

        # Calculate dynamic width based on number of samples for better spacing
        num_samples = len(total_copy_number)
        min_width = 650
        width_per_sample = 20  # Width per sample to ensure good spacing
        calculated_width = max(min_width, num_samples * width_per_sample)
        max_width = 3000  # Maximum width to prevent excessively wide charts
        chart_width = min(calculated_width, max_width)

        # Create bar chart
        bar_chart = alt.Chart(total_copy_number).mark_bar(
            cornerRadius=12,  
            size=22           
        ).encode(
            x=alt.X(
                'Sample', 
                sort=x_sort, 
                axis=alt.Axis(
                    labelFontWeight='bold', 
                    labelColor='#4B5563',
                    labelFontSize=20,  
                    titleFontWeight='bold', 
                    titleColor='#374151',
                    titleFontSize=32,  
                    labelAngle=90,
                    labelOverlap=False, 
                    labelLimit=0,
                    labelPadding=10  
                ),
                # Increase spacing between ticks and bars
                scale=alt.Scale(paddingInner=0.9, paddingOuter=0.2)  
            ),
            y=alt.Y('Total Copy Number', axis=alt.Axis(
                labelFontWeight='bold',
                labelColor='#4B5563',
                labelFontSize=20,   
                titleFontWeight='bold',
                titleColor='#374151',
                titleFontSize=32   
            )),
            tooltip=['Sample', 'Total Copy Number'],
            color=alt.Color('Sample', scale=alt.Scale(
                domain=list(color_mapping.keys()), 
                range=list(color_mapping.values())
            ), legend=None)
        ).properties(
            width=chart_width,    
            height=560,   
            title=alt.TitleParams(
                text='Total Copy Number per Sample',
                fontSize=28,
                fontWeight='bold',
                anchor='middle',
                color='#1F2937',
                subtitle=subtitle_text if subtitle_text else None,
                subtitleFontSize=20,
                subtitleColor='#4c1d95',
                subtitleFontWeight='bold',
                subtitlePadding=25,
            )
        )

        if pathogenic_threshold > 0:
            threshold_line = alt.Chart(pd.DataFrame({'Total Copy Number': [pathogenic_threshold]})).mark_rule(
                color='#EF4444', strokeDash=[5, 5], size=3
            ).encode(y='Total Copy Number:Q')

            threshold_pointer = alt.Chart(pd.DataFrame({'Total Copy Number': [pathogenic_threshold]})).mark_text(
                text='🚨 Pathogenic Threshold', align='left', dx=7, dy=-16, fontSize=24, color='#DC2626', fontWeight='bold'
            ).encode(y='Total Copy Number:Q')
            
            bar_chart = alt.layer(bar_chart, threshold_line, threshold_pointer).configure_view(
                strokeWidth=0,
                fill='rgba(255,255,255,0.9)'
            )
        else:
            bar_chart = bar_chart.configure_view(
                strokeWidth=0,
                fill='rgba(255,255,255,0.9)'
            )

        st.altair_chart(bar_chart, use_container_width=True)

        #add a gap for plots seperation 
        st.markdown("<div style='height: 40px;'></div>", unsafe_allow_html=True)


def stack_plot(record, motif_names, sequences, span_list, motif_ids_list, sort_by="Value", max_height=None, max_width=None):

        # Create stack plot
        motif_colors = utils.get_color_palette(len(record['motifs']))


        motif_colors = {idx: color for idx, color in enumerate(motif_colors)}

        updated_region = record['id']
        df = parsers.create_motif_dataframe(sequences, motif_colors, motif_ids_list, span_list, motif_names, parsers.parse_motif_range)
        if df.empty:
            return motif_colors, df
        df['Length'] = df['End'] - df['Start']

        default_height = 1000 
        chart_height = max(default_height, len(sequences) * 20)
        st.write(chart_height)

        # Apply max_height limit
        if max_height is not None:
            chart_height = min(chart_height, max_height)

        df['Sample'] = df['Sample'].apply(lambda x: x.replace("_pathogenic", ""))
        
        # Calculate total length per sample
        sample_lengths = df.groupby('Sample')['Length'].sum().sort_values(ascending=False)
        
        # Apply sorting
        if sort_by == "Value":
            sorted_samples = sample_lengths.index.tolist()
        else:
            # Sort alphabetically
            sorted_samples = sorted(df['Sample'].unique())
        
        # Convert to categorical
        df['Sample'] = pd.Categorical(df['Sample'], categories=sorted_samples, ordered=True)
        df = df.sort_values(by=['Sample', 'Start', 'End'])
        df['Order'] = range(len(df))

        # Filter out interruptions
        df_filtered = df[df['Motif'] != 'Interruption']

        # Calculate total copy number per sample
        for sample in df_filtered['Sample'].unique():
            df.loc[df['Sample'] == sample, 'Total Copy Number'] = df[df['Sample'] == sample].shape[0]

        # Calculate min and max copy number
        min_copy_number = df['Total Copy Number'].min()
        max_copy_number = df['Total Copy Number'].max()

        # Extract pathogenic threshold information
        pathogenic_threshold = 0
        gene_name = None
        inheritance = None
        disease = None
        above_threshold_samples = None

        # Check if region is pathogenic
        if updated_region in st.session_state.pathogenic_TRs['region'].unique():
            pathogenic_threshold = st.session_state.pathogenic_TRs.loc[
                (st.session_state.pathogenic_TRs['region'] == updated_region)
            ]['pathogenic_min'].values[0]
            if pathogenic_threshold is None:
                pathogenic_threshold = 0
            gene_name = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['gene'].values[0]
            inheritance = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['inheritance'].values[0]
            disease = st.session_state.pathogenic_TRs.loc[
                st.session_state.pathogenic_TRs['region'] == updated_region
            ]['disease'].values[0]

        title = st.empty()

        # Display min and max copy number
        title.markdown(f"""
            <div style="display: flex; gap: 20px; align-items: center; margin-bottom: 15px;">
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="font-size: 20px; color: #64748b; font-weight: 600;">Min:</div>
                    <div style="font-size: 22px; color: #1e293b; font-weight: 700; background: #f1f5f9; padding: 6px 14px; border-radius: 8px; border: 1px solid #e2e8f0;">{min_copy_number}</div>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="font-size: 20px; color: #64748b; font-weight: 600;">Max:</div>
                    <div style="font-size: 22px; color: #1e293b; font-weight: 700; background: #f1f5f9; padding: 6px 14px; border-radius: 8px; border: 1px solid #e2e8f0;">{max_copy_number}</div>
                </div>
            </div>
        """, unsafe_allow_html=True)

        # Create heatmap data
        heatmap_data = df[df['Motif'] != 'Interruption'].groupby(['Sample', 'Motif']).size().reset_index(name='Count')
        
        # Convert to categorical
        heatmap_data['Sample'] = pd.Categorical(heatmap_data['Sample'], categories=sorted_samples, ordered=True)

        # Sort heatmap data by sample
        heatmap_data = heatmap_data.sort_values('Sample')

        # Define the sorting for the heatmap and stack plot
        y_sort = sorted_samples  

        # Prepare stack plot data
        df["pathogenic"] = df["Total Copy Number"].apply(lambda x: "Pathogenic" if x >= pathogenic_threshold else "Not Pathogenic")
        
        # Calculate the pathogenic threshold length
        pathogenic_threshold_length = pathogenic_threshold * len(motif_names[0])

        # Calculate the sequence length per sample
        df['Sequence_length'] = df.groupby('Sample')['Length'].transform('sum')

        # Create heatmap
        motif_count = len(heatmap_data['Motif'].unique())
        width_per_motif = 35
        min_width = 40
        max_width = 900

        # Calculate the dynamic width of the heatmap
        dynamic_width = max(min_width, min(width_per_motif * motif_count, max_width))

        # Create heatmap
        heatmap = alt.Chart(heatmap_data).mark_rect(
            cornerRadius=4,
            stroke='white',
            strokeWidth=1
        ).encode(
            y=alt.Y('Sample:N', 
                title='', 
                sort=y_sort,  
                axis=alt.Axis(
                    labelFontWeight='bold', 
                    labelColor='#4B5563',
                    labelFontSize=28,
                    titleFontWeight='bold',
                    titleColor='#374151',
                    labelLimit=0,
                    ticks=False,
                    tickSize=20,
                    tickOffset=-2,
                    labelPadding=10,
                    labelOverlap=False,
                ),
                scale=alt.Scale(paddingInner=0, paddingOuter=0.6)),
            x=alt.X('Motif:N', 
                title='', 
                sort=alt.EncodingSortField(field='Count', op='sum', order='descending'),
                axis=alt.Axis(
                    labelFontWeight='bold',
                    labelFontSize=28,
                    labelColor='#4B5563',
                    titleFontWeight='bold',
                    titleColor='#374151',
                    labelAngle=90,
                    labelPadding=10,
                    labelLimit=100
                )),
            color=alt.Color('Count:Q', 
                        scale=alt.Scale(scheme='redpurple'), 
                        title='',
                        legend=alt.Legend(
                            title="",
                            orient="top",
                            titleFontSize=28,
                            titleColor='#1F2937',
                            gradientLength=120,
                            labelFontSize=18,
                            labelColor='#4B5563',
                            columns=4,
                            offset=-0,  # decrease offset to bring legend nearer
                            padding=10  # reduce outer padding
                        )),
            tooltip=['Sample', 'Motif', 'Count']
        ).properties(
            width=dynamic_width,
            height=chart_height,
            title=''
        )

        # Create stack plot
        stack_chart = alt.Chart(df).mark_bar(
            cornerRadius=0
        ).encode(
            y=alt.Y(
                'Sample:N',  
                sort=y_sort,  
                title='',
                axis=alt.Axis(
                    labelOverlap=False, 
                    ticks=False, 
                    labelColor='#4B5563', 
                    labelFontSize=0, 
                    labelFontWeight='bold',
                    titleColor='#374151',
                    titleFontSize=36,
                    labelPadding = 1000,
                    labels = False,
                    labelLimit = 0,
                ),
                scale=alt.Scale(paddingInner=0, paddingOuter=0.6)
            ),
            x=alt.X('Length', 
                title='Sequence Length', 
                stack='zero', 
                axis=alt.Axis(
                    labelColor='#4B5563', 
                    labelFontSize=28, 
                    labelFontWeight='bold', 
                    titleColor='#374151',
                    titleFontSize=36,
                    labelOverlap=False,
                    labelAngle=90,  
                    tickMinStep=100,
                ), 
                scale=alt.Scale(domain=[0, df['Sequence_length'].max() + 10])),
            color=alt.Color(
                'Motif', 
                scale=alt.Scale(
                    domain=list(motif_names) + ['Interruption'], 
                    range=list(motif_colors.values()) + ['#EF4444']
                ),
                legend=alt.Legend(
                    title="",
                    orient="top",
                    gradientLength=120,
                    labelFontSize=28,
                    labelColor='#4B5563',
                    symbolStrokeWidth=12,
                    symbolSize=280,
                    symbolType="square",
                    columns=len(motif_names) + 1 if len(motif_names) + 1 <= 4 else 3,
                    columnPadding=15
                )
            ),
            order=alt.Order('Order', sort='ascending'),
            tooltip=['Sample', 'Motif', 'Start', 'End', 'Sequence', 'pathogenic', 'Length', 'Sequence_length']
        ).properties(
            width=max_width if max_width is not None else 'container',
            height=chart_height,
            title='',
        )

        # Create pathogenic threshold elements
        threshold_elements = []
        
        # Calculate above_threshold_samples
        above_threshold_samples = df.groupby('Sample')['Length'].sum().reset_index()
        if pathogenic_threshold > 0:
            above_threshold_samples = above_threshold_samples[above_threshold_samples['Length'] > pathogenic_threshold_length]
        else:
            above_threshold_samples = above_threshold_samples.iloc[0:0]  

        if pathogenic_threshold > 0 and df.groupby('Sample')['Length'].sum().max() > pathogenic_threshold_length:
            # Add threshold line to the stack plot
            rule = alt.Chart(pd.DataFrame({'x': [pathogenic_threshold_length]})).mark_rule(
                color='#EF4444', strokeDash=[5, 5], strokeWidth=2
            ).encode(x='x:Q')
            threshold_elements.append(rule)
            
            if not above_threshold_samples.empty:
                # Create arrow data
                arrow_data = above_threshold_samples.copy()
                arrow_data['arrow_x'] = arrow_data['Length'] + (df['Sequence_length'].max() * 0.01)  
                arrow_data['arrow_y_offset'] = 0.4  
                

                arrows = alt.Chart(arrow_data).mark_text(
                    text='➔',  
                    align='left',
                    baseline='middle',
                    fontSize=28,
                    color='#DC2626',
                    fontWeight='bold',
                    angle=235,   
                    dy=0,
                    dx=-10        
                ).encode(
                    x=alt.X('arrow_x:Q', title=None),
                    y=alt.Y('Sample:N', sort=y_sort),
                    tooltip=[alt.Tooltip('Sample:N'), alt.Tooltip('Length:Q', title='Total Length')]
                )
                
                threshold_elements.append(arrows)
                
                # Add warning text for pathogenic samples to the stack plot
                warning_text = alt.Chart(arrow_data).mark_text(
                    text=' PATHOGENIC',
                    align='left',
                    baseline='middle',
                    fontSize=18,
                    color='#DC2626',
                    fontWeight='bold',
                    dy=10,
                    dx=10,
                ).encode(
                    x=alt.X('arrow_x:Q', title=None),
                    y=alt.Y('Sample:N', sort=y_sort),
                    tooltip=[alt.Tooltip('Sample:N'), alt.Tooltip('Length:Q', title='Total Length')]
                )
                
                threshold_elements.append(warning_text)

        # Combine all elements
        if threshold_elements:
            all_stack_elements = [stack_chart] + threshold_elements
            final_stack_chart = alt.layer(*all_stack_elements)
        else:
            final_stack_chart = stack_chart

        # Create chart title
        chart_title = "Motif Distribution across Samples"
        subtitle_text = ""
        if gene_name or inheritance or disease or updated_region:
            gene_info_parts = []
            if gene_name:
                gene_info_parts.append(f"Gene: {gene_name}")
            if inheritance:
                gene_info_parts.append(f"Inheritance: {inheritance}")
            if disease:
                disease_display = disease 
                gene_info_parts.append(f"Disease: {disease_display}")
            if updated_region:
                gene_info_parts.append(f"Region: {updated_region}")
            subtitle_text = " • ".join(gene_info_parts)

        min_heatmap_display_width = 120  

      
        if dynamic_width == 40:
            # Pad the heatmap with a blank chart to the left
            padding_width = min_heatmap_display_width - dynamic_width

            blank_left = alt.Chart(
                pd.DataFrame({"dummy": [0]})
            ).mark_rect(opacity=0).encode(
                x=alt.value(0),
                y=alt.value(0)
            ).properties(
                width=padding_width,
                height=chart_height
            )

            heatmap_display = alt.hconcat(
                blank_left,
                heatmap.properties(width=dynamic_width),
                spacing=0
            )
        else:
            # Display the heatmap
            heatmap_display = heatmap

        combined_chart = alt.hconcat(
            heatmap_display, 
            final_stack_chart,
            spacing=0  
        ).resolve_scale(
            y='shared'
        ).properties(
            title=alt.TitleParams(
                text=chart_title,
                fontSize=32,
                fontWeight='bold',
                anchor='middle',
                color='#1F2937',
                subtitle=subtitle_text,
                subtitleFontSize=20,
                subtitleColor='#4c1d95',
                subtitleFontWeight='bold',
                subtitlePadding=25,
            )
        ).configure_view(
            strokeWidth=0
        ).configure_scale(
            bandPaddingInner=0
        ).configure_axis(
            labelLimit=1000,
            tickCount=len(df['Sample'].unique()),
            labelOverlap=False,
            tickMinStep=max(200, int(df['Sequence_length'].max() * 0.1)),
        ).configure_axisY(
            labelFontSize=12,
            labelLimit=1000,
            tickCount=len(df['Sample'].unique()),
            labelOverlap=False,            
        )

        
        if chart_height > default_height:
            combined_chart = combined_chart.configure_axisY(labelFontSize=8)
        
        # Display the combined chart
        st.altair_chart(combined_chart, use_container_width=True)
        
        # Add spacing after stack plot to prevent overlap with next plot title
        st.markdown("<div style='height: 40px;'></div>", unsafe_allow_html=True)
        
        # Show pathogenic info if applicable
        if pathogenic_threshold > 0 and above_threshold_samples is not None and not above_threshold_samples.empty:
            st.warning(f"🚨 **Pathogenic Alert**: {len(above_threshold_samples)} sample(s) exceed the pathogenic threshold of {pathogenic_threshold} copies")
        
        return motif_colors, df

def plot_Cohort_results(cohort_records, cohort_mode):
        
        sequences = []
        span_list = []
        motif_ids_list = []
        # make space between the last print 
        sort_by = st.radio("Sort by:", ("Value", "Sample Name"), horizontal=True, key="sort_by_cohort")
        
        # Extract genotype information for comparison
        genotypes_dict = {}


        if cohort_mode == "assembly":
            # For assembly mode, group samples by base name and combine haplotypes
            sample_groups = {}
            for key in cohort_records.keys():
                # Check if sample name ends with _h1 or _h2
                if key.endswith('_h1'):
                    base_name = key[:-3]  # Remove _h1 suffix
                    if base_name not in sample_groups:
                        sample_groups[base_name] = {}
                    sample_groups[base_name]['h1'] = cohort_records[key]
                elif key.endswith('_h2'):
                    base_name = key[:-3]  # Remove _h2 suffix
                    if base_name not in sample_groups:
                        sample_groups[base_name] = {}
                    sample_groups[base_name]['h2'] = cohort_records[key]
                else:
                    # Sample doesn't have _h1/_h2 suffix, use as is
                    genotypes_dict[key] = cohort_records[key]['gt']
            
            # Compute diploid genotypes for grouped samples
            for base_name, haplotypes in sample_groups.items():
                if 'h1' in haplotypes and 'h2' in haplotypes:
                    h1_record = haplotypes['h1']
                    h2_record = haplotypes['h2']
                    
                    # Extract genotypes and IDs
                    gt_h1 = h1_record['gt']
                    gt_h2 = h2_record['gt']
                    ids_h1 = h1_record.get('motif_ids_h', [])
                    ids_h2 = h2_record.get('motif_ids_h', [])
                    
                    # Compute diploid genotype
                    diploid_gt = parsers.compute_diploid_genotype_assembly(gt_h1, gt_h2, ids_h1, ids_h2)
                    genotypes_dict[base_name] = diploid_gt
                elif 'h1' in haplotypes:
                    # Only h1 available
                    genotypes_dict[base_name] = haplotypes['h1']['gt']
                elif 'h2' in haplotypes:
                    # Only h2 available
                    genotypes_dict[base_name] = haplotypes['h2']['gt']
        else:
            # For reads mode, use genotypes as is
            for key in cohort_records.keys():
                genotypes_dict[key] = cohort_records[key]['gt']
        
        # Display genotype comparison matrix
        st.markdown("---")
        create_genotype_comparison_matrix(genotypes_dict)
        st.markdown("---")
        if cohort_mode == "assembly":
            for key in cohort_records.keys():
                sequences.append({'name': f'{key}', 'sequence': cohort_records[key]['alt_allele']})
                span_list.append(cohort_records[key]['spans'])
                motif_ids_list.append(cohort_records[key]['motif_ids_h'])
        else:
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

        motif_colors, df = stack_plot(record, motif_names, sequences, span_list, motif_ids_list, sort_by)

        region = f"{record['chr']}:{record['pos']-1}-{record['stop']-1}"
        bar_plot_motif_count(df, region, sort_by=sort_by)




        figure = go.Figure()

        if df.empty:
            st.warning("No motifs found in the region")
            st.stop()
        unique_samples = df['Sample'].unique()
        unique_samples = [sample for sample in unique_samples if sample != "Interruption"]

        ref_data = []
        allele_data = []

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
                    marker=dict(size=24)  
                )
                if sample in ["Ref", "Allel1", "Allel2"]:
                    if sample == "Ref":
                        trace['marker']['symbol'] = "x"
                        ref_data.append(trace)
                    else:
                        trace['marker']['symbol'] = "triangle-down"
                        allele_data.append(trace)
                else:
                    figure.add_trace(trace)

        for trace in ref_data + allele_data:
            figure.add_trace(trace)
        color_mapping = { sample: np.random.choice(px.colors.qualitative.Plotly) for sample in unique_samples}
        color_mapping["Ref"] = "green"

        figure.for_each_trace(lambda trace: trace.update(marker=dict(color=color_mapping.get(trace.name, 'gray'))))

        unique_legend_names = set()
        for trace in figure.data:
            if trace.name in unique_legend_names:
                trace.showlegend = False
            elif trace.name in color_mapping:
                unique_legend_names.add(trace.name)
                trace.showlegend = True
            else:
                trace.showlegend = False

        big_font_size = 36      
        axis_font_size = 22
        legend_font_size = 22
        tick_font_size = 20
        font_family = 'Arial'
        font_weight = 'bold'

        figure.update_layout(
            title=dict(
                text="<b>Motif Occurrences</b>",
                font=dict(size=big_font_size, family=font_family, color='black',),
                x=0.5,
                xanchor='center'
            ),
            xaxis_title=dict(
                text="<b>Motif</b>",
                font=dict(size=axis_font_size, family=font_family, color='black')
            ),
            yaxis_title=dict(
                text="<b>Count</b>",
                font=dict(size=axis_font_size, family=font_family, color='black')
            ),
            font=dict(
                size=tick_font_size,
                family=font_family,
                color='black'
            ),
            legend=dict(
                font=dict(
                    size=legend_font_size,
                    family=font_family,
                    color='black'
                )
            )
        )

        xaxis_colors = {motif: motif_colors[idx] for idx, motif in enumerate(motif_names)}
        figure.update_xaxes(
            tickmode='array',
            tickvals=list(xaxis_colors.keys()),
            ticktext=[
                f'<span style="color:{xaxis_colors[motif]}; font-size:{tick_font_size+4}px; font-weight:bold">{motif}</span>'
                for motif in xaxis_colors.keys()
            ],
            tickangle=45,
            tickfont=dict(size=tick_font_size, family=font_family, color="black",),
        )
        figure.update_yaxes(
            range=[0, df['Sample'].value_counts().max()],
            tickfont=dict(size=tick_font_size, family=font_family, color="black",),
        )

        st.plotly_chart(figure, use_container_width=True)



def plot_HGSVC_VS_allele(record, hgsvc_records, motif_names):
        
        """
        This function plots the HGSVC vs allele plots
        """

        sequences = []
        span_list = []
        motif_ids_list = []
        
        
        # Radio group display for sorting the plots
        st.html("""
            <style>
                /* Outer radio styling */
                .custom-radio-group {
                    display: flex;
                    gap: 22px;
                    align-items: center;
                    background: rgba(255, 255, 255, 0.97);
                    border-radius: 16px;
                    padding: 18px 32px 18px 20px;
                    border: 1.5px solid #c7d2fe;
                    box-shadow: 0 2px 16px rgba(102, 126, 234, 0.07);
                    margin: 22px 0 18px 0;
                    backdrop-filter: blur(14px);
                }
                .custom-radio-label {
                    display: flex !important;
                    align-items: center;
                    gap: 7px;
                    margin-right: 0 !important;
                }
                .custom-radio-btn {
                    appearance: none;
                    width: 22px;
                    height: 22px;
                    border-radius: 50%;
                    border: 2.4px solid #bfc8fa;
                    transition: border .2s, box-shadow .2s;
                    outline: none;
                    margin-right: 12px;
                    box-sizing: border-box;
                    background: white;
                    box-shadow: 0 0 0px 4px #e1eaff00;
                    position: relative;
                }
                .custom-radio-btn:checked {
                    border-color: #667eea;
                    background: linear-gradient(135deg, #667eea 60%, #764ba2 100%);
                    box-shadow: 0 0 0px 4px #bfc8fa44;
                }
                .custom-radio-btn:checked::before {
                    content: '';
                    position: absolute;
                    left: 6.5px; top: 6.5px;
                    width: 7px; height: 7px;
                    border-radius: 50%;
                    background: #fff;
                }
                .custom-radio-label.active,
                .custom-radio-label:hover {
                    background: linear-gradient(135deg, #f3f5fa 90%, #eef0fe 100%);
                    border-radius: 9px;
                    box-shadow: 0 2px 12px #e1eaff33;
                    cursor: pointer;
                }
                .custom-radio-label .radio-emoji {
                    font-size: 1.16em;
                }
                .sortby-title {
                    font-size: 19px;
                    font-weight: 700;
                    color: #43486b;
                    margin-right: 9px;
                }
            </style>
        """)

        sortby_options = [
            {"label": '<span class="radio-emoji">📊</span> Value', "value": "📊 Value"},
            {"label": '<span class="radio-emoji">👥</span> Sample Name', "value": "👥 Sample Name"}
        ]

        # Maintain radio state with Streamlit session_state
        if "sort_by_cohort" not in st.session_state:
            st.session_state["sort_by_cohort"] = sortby_options[0]["value"]

        st.html(
            '<div class="custom-radio-group">' +
            '<span class="sortby-title">Sort by:</span>' +
            ''.join(
                f'''<label class="custom-radio-label{' active' if st.session_state["sort_by_cohort"] == option['value'] else ''}">
                    <input type="radio" class="custom-radio-btn" name="sort_by_cohort" value="{option['value']}" 
                        {'checked' if st.session_state["sort_by_cohort"] == option['value'] else ''} 
                        onclick="window.dispatchEvent(new CustomEvent('sortbyChange', {{detail: '{option['value']}'}}))">
                    {option['label']}
                </label>''' 
                for option in sortby_options
            ) +
            '</div>',
            
        )

        # Custom events solution to handle client-side radio change
        st.components.v1.html(f"""
        <script>
            window.addEventListener('sortbyChange', function(e) {{
                // Use a cookie to communicate - Streamlit will reload when widget is redrawn
                document.cookie = "sortby_cohort_choice=" + encodeURIComponent(e.detail) + "; path=/";
                window.location.reload();
            }});
        </script>
        """, height=0)

        motif_colors = utils.get_color_palette(len(motif_names))

        # Server-side get cookie function
        def get_cookie_value(key):
            try:
                cookies_param = st.query_params().get('__streamlit_cookies__')
                if cookies_param:
                    # Handle both string and list cases
                    if isinstance(cookies_param, list):
                        cookies = cookies_param[0] if cookies_param else ''
                    else:
                        cookies = str(cookies_param)
                else:
                    cookies = ''
                m = _re.search(rf'{_re.escape(key)}=(.*?)(?:;|$)', cookies)
                return m.group(1) if m else None
            except Exception:
                return None

        # Streamlit doesn't provide cookie parsing directly, so we try with components
        sortby_cookie = get_cookie_value('sortby_HGSVC_choice')
        if sortby_cookie in [o["value"] for o in sortby_options]:
            st.session_state["sort_by_HGSVC"] = sortby_cookie

        sort_by = st.session_state.get("sort_by_HGSVC", sortby_options[0]["value"])
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
            span_list.append(hgsvc_records[key]['spans'])
            motif_ids_list.append(hgsvc_records[key]['motif_ids_h'])
        motif_colors, df = stack_plot(record, motif_names, sequences, span_list, motif_ids_list, sort_by, None, None)
        region = f"{record['chr']}:{record['pos']-1}-{record['stop']-1}"
        
        # Add spacing before bar plot to prevent overlap
        st.markdown("<div style='height: 10px;'></div>", unsafe_allow_html=True)
        
        # Create two columns for side-by-side plots
        st.markdown('<div class="plot-card">', unsafe_allow_html=True)
        bar_plot_motif_count(df, region, sort_by)
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Add spacing before cluster plot to prevent overlap
        st.markdown("<div style='height: 40px;'></div>", unsafe_allow_html=True)
        
        # Cluster plot in full width below
        st.markdown('<div class="plot-card-full">', unsafe_allow_html=True)
        cluster_plot_samples(df, motif_colors)
        st.markdown('</div>', unsafe_allow_html=True)


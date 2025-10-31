"""
Alternative visualization methods for genotype display instead of st.expander
"""

import streamlit as st
import plotly.graph_objects as go
import pandas as pd

def create_genotype_with_tabs(genotypes_dict):
    """Using Streamlit tabs for genotype display"""
    
    samples = list(genotypes_dict.keys())
    unique_genotypes = list(set(genotypes_dict.values()))
    
    # Group samples by genotype
    genotype_groups = {}
    for sample, gt in genotypes_dict.items():
        if gt not in genotype_groups:
            genotype_groups[gt] = []
        genotype_groups[gt].append(sample)
    
    # Create tabs for each unique genotype
    tabs = st.tabs([f"{gt} ({len(genotype_groups[gt])})" for gt in unique_genotypes])
    
    for idx, (gt, samples_list) in enumerate(genotype_groups.items()):
        with tabs[idx]:
            interpretation = interpret_genotype(gt)
            
            # Header with color coding
            st.markdown(f"""
                <div style="background: linear-gradient(135deg, {interpretation['color']}22, {interpretation['bg_color']});
                            padding: 20px; border-radius: 10px; margin-bottom: 20px;
                            border-left: 5px solid {interpretation['color']};">
                    <h2>{interpretation['icon']} {interpretation['status'].replace('_', ' ').title()}</h2>
                    <p style="font-size: 24px; font-weight: bold; color: {interpretation['color']};">
                        Genotype: {gt}
                    </p>
                </div>
            """, unsafe_allow_html=True)
            
            # Show samples in this group
            st.markdown(f"### Samples: {len(samples_list)}")
            for sample in samples_list:
                st.markdown(f"- **{sample}**")
            
            # Create a mini-chart
            if len(samples_list) > 1:
                fig = go.Figure()
                fig.add_trace(go.Bar(
                    x=[sample for sample in samples_list],
                    y=[1] * len(samples_list),
                    marker_color=interpretation['color']
                ))
                fig.update_layout(height=200, showlegend=False, 
                                title=f"Samples with genotype {gt}")
                st.plotly_chart(fig, use_container_width=True)


def create_genotype_with_columns(genotypes_dict):
    """Using responsive columns for genotype display"""
    
    samples = list(genotypes_dict.keys())
    unique_genotypes = list(set(genotypes_dict.values()))
    
    st.markdown("### 🧬 Genotype Distribution")
    
    # Create columns based on number of unique genotypes
    num_cols = min(len(unique_genotypes), 4)  # Max 4 columns
    cols = st.columns(num_cols)
    
    for idx, gt in enumerate(unique_genotypes):
        col = cols[idx % num_cols]
        
        with col:
            interpretation = interpret_genotype(gt)
            samples_with_gt = [s for s, g in genotypes_dict.items() if g == gt]
            
            st.markdown(f"""
                <div style="background: {interpretation['bg_color']};
                            border: 2px solid {interpretation['color']};
                            border-radius: 10px; padding: 15px; margin-bottom: 10px;">
                    <h4 style="color: {interpretation['color']};">{interpretation['icon']}</h4>
                    <h3 style="margin: 0;">{gt}</h3>
                    <p style="color: #666; font-size: 0.9em;">{len(samples_with_gt)} samples</p>
                </div>
            """, unsafe_allow_html=True)


def create_genotype_with_popover(genotypes_dict):
    """Using a compact view with modal-like display"""
    
    st.markdown("### 🧬 Genotype Overview")
    
    # Summary stats
    unique_genotypes = list(set(genotypes_dict.values()))
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Samples", len(genotypes_dict))
    with col2:
        st.metric("Unique Genotypes", len(unique_genotypes))
    with col3:
        st.metric("Most Common", max(set(genotypes_dict.values()), 
                                     key=genotypes_dict.values().count))
    
    # Compact genotype cards
    if st.button("📋 View Detailed Genotypes", use_container_width=True):
        st.markdown("---")
        st.markdown("### Detailed Genotype Information")
        
        # Group by genotype
        genotype_groups = {}
        for sample, gt in genotypes_dict.items():
            if gt not in genotype_groups:
                genotype_groups[gt] = []
            genotype_groups[gt].append(sample)
        
        for gt, samples_list in genotype_groups.items():
            interpretation = interpret_genotype(gt)
            
            with st.container():
                st.markdown(f"""
                    <div style="background: {interpretation['bg_color']};
                                border-left: 4px solid {interpretation['color']};
                                padding: 15px; border-radius: 5px; margin: 10px 0;">
                """, unsafe_allow_html=True)
                
                st.markdown(f"**{interpretation['icon']} {gt}** - *{interpretation['status'].replace('_', ' ').title()}*")
                
                # Show samples in this group
                sample_text = ", ".join(samples_list)
                st.markdown(f"*Samples: {sample_text}*")
                
                st.markdown("</div>", unsafe_allow_html=True)


def create_genotype_with_chart(genotypes_dict):
    """Using an interactive chart for genotype visualization"""
    
    # Count genotypes
    genotype_counts = {}
    for gt in genotypes_dict.values():
        genotype_counts[gt] = genotype_counts.get(gt, 0) + 1
    
    # Create bar chart
    fig = go.Figure()
    
    colors = []
    labels = []
    for gt, count in genotype_counts.items():
        interpretation = interpret_genotype(gt)
        colors.append(interpretation['color'])
        labels.append(f"{gt}<br>{interpretation['status'].replace('_', ' ').title()}")
    
    fig.add_trace(go.Bar(
        x=list(genotype_counts.keys()),
        y=list(genotype_counts.values()),
        marker_color=colors,
        text=[f"{c} samples" for c in genotype_counts.values()],
        textposition='auto',
        hovertemplate='<b>%{x}</b><br>%{y} samples<extra></extra>'
    ))
    
    fig.update_layout(
        title="Genotype Distribution",
        xaxis_title="Genotype",
        yaxis_title="Number of Samples",
        height=400,
        hovermode='closest'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Show details in an expandable section
    with st.expander("📋 View Sample Details", expanded=False):
        for sample, gt in sorted(genotypes_dict.items()):
            interpretation = interpret_genotype(gt)
            st.markdown(f"**{sample}**: {gt} - {interpretation['status'].replace('_', ' ').title()}")


def create_genotype_with_table(genotypes_dict):
    """Using an interactive data table for genotype display"""
    
    st.markdown("### 🧬 Genotype Table")
    
    # Create DataFrame
    df = pd.DataFrame([
        {'Sample': sample, 'Genotype': gt} 
        for sample, gt in genotypes_dict.items()
    ])
    
    # Add interpretation
    df['Status'] = df['Genotype'].apply(
        lambda x: interpret_genotype(x)['status'].replace('_', ' ').title()
    )
    df['Color'] = df['Genotype'].apply(
        lambda x: interpret_genotype(x)['color']
    )
    
    # Use st.dataframe with formatting
    st.dataframe(
        df[['Sample', 'Genotype', 'Status']],
        use_container_width=True,
        hide_index=True
    )
    
    # Group summary
    st.markdown("### Summary by Genotype")
    summary = df.groupby('Genotype').size().reset_index(name='Count')
    st.dataframe(summary, use_container_width=True, hide_index=True)


def create_genotype_with_slider(genotypes_dict):
    """Interactive genotype viewer with slider navigation"""
    
    samples = list(genotypes_dict.keys())
    
    st.markdown("### 🧬 Sample Genotype Viewer")
    
    # Slider to navigate through samples
    sample_idx = st.slider("Select Sample", 0, len(samples)-1, 0)
    
    selected_sample = samples[sample_idx]
    selected_genotype = genotypes_dict[selected_sample]
    interpretation = interpret_genotype(selected_genotype)
    
    # Display current sample
    st.markdown(f"""
        <div style="background: {interpretation['bg_color']};
                    border: 3px solid {interpretation['color']};
                    border-radius: 15px; padding: 30px; text-align: center;">
            <h2>{interpretation['icon']}</h2>
            <h1 style="color: {interpretation['color']};">{selected_sample}</h1>
            <h3>Genotype: <strong>{selected_genotype}</strong></h3>
            <p style="font-size: 1.2em;">{interpretation['status'].replace('_', ' ').title()}</p>
        </div>
    """, unsafe_allow_html=True)
    
    # Navigation info
    col1, col2 = st.columns(2)
    with col1:
        if st.button("⬅️ Previous"):
            st.session_state.sample_idx = (sample_idx - 1) % len(samples)
    with col2:
        if st.button("Next ➡️"):
            st.session_state.sample_idx = (sample_idx + 1) % len(samples)


# Helper function (you'll need to implement this based on your existing code)
def interpret_genotype(gt):
    """Interpret genotype and return color/status/icon"""
    # Placeholder - implement based on your actual logic
    if gt.startswith('Homozygous'):
        return {'color': '#10B981', 'bg_color': '#D1FAE5', 
                'status': 'homozygous', 'icon': '⭕'}
    elif gt.startswith('Heterozygous'):
        return {'color': '#3B82F6', 'bg_color': '#DBEAFE', 
                'status': 'heterozygous', 'icon': '🔵'}
    else:
        return {'color': '#F59E0B', 'bg_color': '#FEF3C7', 
                'status': 'unknown', 'icon': '❓'}


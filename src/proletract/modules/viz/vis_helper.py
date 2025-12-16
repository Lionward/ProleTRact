import pandas as pd 
import re
import streamlit as st
import altair as alt
import hashlib
import html
from collections import Counter

def apply_global_typography() -> None:
    key = "pt_global_typography_v4"
    if st.session_state.get(key):
        return

    st.markdown(
        """
        <style>
        :root {
            --pt-font-base: 2px;
            --pt-font-xs: 0.78rem;
            --pt-font-sm: 0.9rem;
            --pt-font-md: 1rem;
            --pt-font-lg: 0.15rem;
            --pt-font-xl: 0.32rem;
            --pt-font-xxl: 0.55rem;
            --pt-weight-regular: 400;
            --pt-weight-semibold: 600;
            --pt-weight-bold: 700;
        }
        body, .stApp {
            font-family: "Inter", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
            font-size: var(--pt-font-base);
            line-height: 1.58;
            color: #1f2937;
        }
        h1 { font-size: var(--pt-font-xxl); font-weight: var(--pt-weight-bold); }
        h2 { font-size: var(--pt-font-xl); font-weight: var(--pt-weight-semibold); }
        h3 { font-size: var(--pt-font-lg); font-weight: var(--pt-weight-semibold); }
        p, li, label, [data-testid="stWidgetLabel"] > div {
            font-size: var(--pt-font-md);
        }
        .pt-text-xs { font-size: var(--pt-font-xs); }
        .pt-text-sm { font-size: var(--pt-font-sm); }
        .pt-text-md { font-size: var(--pt-font-md); }
        .pt-text-lg { font-size: var(--pt-font-lg); }
        .pt-heading { font-size: var(--pt-font-xl); font-weight: var(--pt-weight-bold); }
        .pt-heading-lg { font-size: var(--pt-font-xxl); font-weight: var(--pt-weight-bold); }
        .stButton > button, .stDownloadButton > button {
            font-size: var(--pt-font-sm);
            font-weight: var(--pt-weight-semibold);
        }
        .stRadio label, .stCheckbox label, .stSelectbox label {
            font-size: var(--pt-font-sm);
        }
        .stMarkdown p, .stMarkdown li {
            font-size: var(--pt-font-md);
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.session_state[key] = True

def parse_motif_range(motif_range):
    pattern = re.compile(r'\((\d+)-(\d+)\)')
    matches = pattern.findall(motif_range)
    ranges = [(int(start)-1, int(end)-1) for start, end in matches]
    return ranges


def motif_legend_html(motif_ids, motif_colors, motif_names):
    """
    Generate a motif legend: improved spacing, font-size 14px, much less vertical stacking,
    but overall much more compact in height.
    """
    apply_global_typography()
    legend_html = ""
    unique_motifs = sorted(set(int(mid) for mid in motif_colors.keys()))
    motif_display_data = []
    for motif_id in unique_motifs:
        color = motif_colors[int(motif_id)]
        motif_name = motif_names[int(motif_id)]
        motif_size = len(motif_name)
        motif_display_data.append({
            'id': motif_id,
            'color': color,
            'name': motif_name,
            'size': motif_size,
        })

    motif_display_data.sort(key=lambda x: x['size'], reverse=True)

    for motif_data in motif_display_data:
        motif_name = motif_data['name']
        display_name = motif_name if len(motif_name) <= 22 else f"{motif_name[:19]}..."
        legend_html += f"""
        <div class="legend-item" data-motif="{motif_data['id']}" title="{motif_name} - {len(motif_name)} bases">
            <span class="legend-color" style="background-color:{motif_data['color']};"></span>
            <span class="legend-text">{display_name} <span class="legend-bps">({len(motif_name)}bp)</span></span>
        </div>
        """

    legend_html += """
        <div class="legend-item" title="Interruption regions between motifs">
            <span class="legend-color interruption-color"></span>
            <span class="legend-text">Interruption</span>
        </div>
    """

    # Restyle for compactness: reduce paddings/margins,
    # single-line motif grid, less header and stats height
    st.markdown("""
        <style>
            .motif-legend-container {
                margin: 0;
                background: #fff;
                border-radius: 8px;
                box-shadow: 0 0.5px 2px rgba(60,80,130,0.06);
                border: 1px solid #e2e8f0;
                width: 100%;
                font-size: 14px !important;
                padding: 7px 12px 7px 12px;
                max-width: 100vw;
            }
            .motif-legend-header {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 4px 10px 3px 10px;
                border-radius: 7px 7px 0 0;
                font-weight: 600;
                font-size: 14px !important;
                display: flex;
                justify-content: space-between;
                align-items: center;
                min-height: 19px;
                margin: -7px -12px 2px -12px;
                line-height: 1.12;
            }
            .motif-count {
                background: rgba(255,255,255,0.19);
                padding: 0 7px;
                border-radius: 6px;
                font-size: 13px !important;
                font-weight: 600;
                height: 18px;
                display: flex;
                align-items: center;
                border: 1px solid #8aa4ed;
            }
            .motif-legend-content {
                padding: 0;
                max-height: 48px;
                overflow-x: auto;
                overflow-y: hidden;
                font-size: 14px !important;
            }
            .motif-legend-grid {
                display: flex;
                flex-wrap: nowrap;
                gap: 8px; /* horizontal gap only */
            }
            .legend-item {
                display: flex;
                align-items: center;
                padding: 1px 6px 1px 3px;
                background: #f8fafc;
                border-radius: 5px;
                border: 1px solid transparent;
                min-height: 22px;
                font-size: 14px !important;
                margin-bottom: 0;
                transition: border 0.1s, background 0.1s;
            }
            .legend-item:hover {
                border-color: #667eea;
                background: #eef1fc;
            }
            .legend-color {
                width: 15px;
                height: 15px;
                border-radius: 50%;
                margin-right: 5px;
                border: 1px solid white;
                flex-shrink: 0;
                box-shadow: 0 1px 2.5px rgba(110,110,150,0.07);
            }
            .interruption-color {
                background: linear-gradient(135deg, #fc8181, #e53e3e) !important;
            }
            .legend-text {
                font-size: 14px !important;
                font-weight: 600 !important;
                color: #374151;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
                max-width: 130px;
                font-family: inherit;
                display: flex;
                align-items: center;
            }
            .legend-bps {
                font-weight: 400;
                font-size: 13px !important;
                color: #667eea;
                margin-left: 5px;
            }
            .legend-stats {
                display: flex;
                gap: 12px;
                padding: 2px 2px 2px 2px;
                background: #eef4ff;
                border-top: 1px solid #e2e8f0;
                font-size: 13px !important;
                color: #485061;
                border-radius: 0 0 7px 7px;
                margin-bottom: 2px;
                margin-top: 1px;
                font-weight: 500;
                min-height: 14px;
            }
            .stat-item {
                display: flex;
                align-items: center;
                gap: 4px;
            }
            .stat-value {
                background: white;
                padding: 0 6px;
                border-radius: 5px;
                border: 1px solid #d1dae8;
                font-weight: 700;
                color: #1f2937;
                font-size: 13px !important;
                height: 18px;
                display: flex;
                align-items: center;
            }
            .motif-legend-content::-webkit-scrollbar {
                height: 6px;
            }
            .motif-legend-content::-webkit-scrollbar-thumb {
                background: #d2dbe7;
                border-radius: 3px;
            }
            .motif-legend-content::-webkit-scrollbar-track {
                background: #f4f6f9;
            }
            @media (max-width: 700px) {
                .legend-text { max-width: 55vw !important; }
            }
        </style>
    """, unsafe_allow_html=True)
    total_motifs = len(unique_motifs)
    motif_sizes = [len(motif_names[int(motif_id)]) for motif_id in unique_motifs]
    avg_size = sum(motif_sizes) / len(motif_sizes) if motif_sizes else 0
    max_size = max(motif_sizes) if motif_sizes else 0
    min_size = min(motif_sizes) if motif_sizes else 0

    st.html(f"""
        <div class="motif-legend-container">
            <div class="motif-legend-header">
                <span>Motifs in Region</span>
                <span class="motif-count">{total_motifs}</span>
            </div>
            <div class="legend-stats">
                <div class="stat-item">
                    <span>Size:</span><span class="stat-value">{min_size}-{max_size}bp</span>
                </div>
                <div class="stat-item">
                    <span>Avg:</span><span class="stat-value">{avg_size:.1f}bp</span>
                </div>
            </div>
            <div class="motif-legend-content">
                <div class="motif-legend-grid">
                    {legend_html}
                </div>
            </div>
        </div>
    """)


def display_dynamic_sequence_with_highlighted_motifs(sequence_name, sequence, motif_ids, spans, motif_colors, motif_names, supporting_reads=None):
    style_config = {
        "fonts": {
            "dashboard": "0.93rem",             
            "header": "0.96rem",               
            "length": "0.91rem",                
            "legend_container": "0.93rem",     
            "legend_header": "0.94rem",        
            "legend_count": "0.91rem",          
            "stats": "0.91rem",                 
        },
        "no_motif": {
            "dashboard_margin": "10px",               
            "container_border_radius": "8px",         
            "header_padding": "5px 8px",              
            "header_font": "10px",                    
            "header_gap": "6px",                      
            "length_font": "10px",                    
            "length_padding": "0px 6px",              
            "length_border_radius": "8px",            
            "length_letter_spacing": "0.4px",         
            "section_padding": "7px",                 
            "section_margin": "7px",                  
            "section_border_radius": "6px",           
            "icon_size": "13px",                      
            "message_font": "13px",                   
            "desc_font": "9px",                       
            "desc_opacity": "0.76",                   
            "content_padding": "6px",                 
            "content_font": "11px",                   
            "content_line_height": "1.4",             
            "scale_height": "30px",                   
            "scale_padding": "7px 7px 10px 7px",      
            "scale_font": "11px",                     
            "scale_min_height": "28px",               
            "scale_marker_offset": "7px",             
        },
        "sequence": {
            "dashboard_margin": "10px",               
            "container_border_radius": "8px",         
            "header_padding": "6px 8px",              
            "header_font": "11px",                    
            "length_font": "11px",                    
            "length_padding": "1px 6px",              
            "length_border_radius": "8px",            
            "stats_gap": "7px",                       
            "stats_padding": "5px 8px",               
            "stats_font": "11px",                     
            "content_padding": "10px",                
            "content_font": "10px",                   
            "content_min_height": "28px",             
            "content_line_height": "1.3",             
            "scale_height": "30px",                   
            "scale_padding": "7px 8px 10px 8px",      
            "scale_font": "11px",                     
            "scale_min_height": "28px",               
            "scale_marker_offset": "11px",            
            "motif_padding": "2px 0px",               
            "motif_border_radius": "7px",             
            "motif_font": "11px",                     
            "interruption_padding": "2px 0px",        
            "interruption_border_radius": "7px",      
        },
        "tooltip": {
            "padding": "7px 10px",                    
            "font": "10px",                           
            "border_radius": "8px",                   
            "max_width": "240px",                     
            "code_padding": "2px 4px",                
            "code_font": "9px",                       
        },
        "scrollbars": {
            "height": "4px"                           
        }
    }
    fonts = style_config["fonts"]
    nm_sizes = style_config["no_motif"]
    seq_sizes = style_config["sequence"]
    tooltip_sizes = style_config["tooltip"]
    scrollbar_sizes = style_config["scrollbars"]
    # make fonts bigger
    st.markdown(f"""
        <style>
            .sequence-dashboard, .sequence-dashboard * {{ font-size: {fonts['dashboard']} !important; }}
            .sequence-header {{ font-size: {fonts['header']} !important; }}
            .sequence-length {{ font-size: {fonts['length']} !important; }}
            .motif-legend-container, .motif-legend-container * {{ font-size: {fonts['legend_container']} !important; }}
            .motif-legend-header {{ font-size: {fonts['legend_header']} !important; }}
            .motif-count {{ font-size: {fonts['legend_count']} !important; }}
            .legend-stats .stat-item span, .legend-stats .stat-item .stat-value {{ font-size: {fonts['stats']} !important; }}
        </style>
    """, unsafe_allow_html=True)

    def make_dom_slug(value: str) -> str:
        slug = re.sub(r'[^A-Za-z0-9_-]+', '-', value).strip('-')
        if not slug:
            slug = "sequence"
        if slug[0].isdigit():
            slug = f"seq-{slug}"
        return re.sub(r'-{2,}', '-', slug)

    dom_slug = make_dom_slug(sequence_name)
    dom_hash = hashlib.md5(sequence_name.encode("utf-8")).hexdigest()[:6]
    sequence_dom_key = f"{dom_slug}-{dom_hash}"
    sequence_content_id = f"sequence-content-{sequence_dom_key}"
    sequence_scale_id = f"sequence-scale-{sequence_dom_key}"
    # Handle cases where motifs are missing or sequence is just one base
    if (motif_ids == ["."]) or (isinstance(sequence, str) and len(sequence) <= 1):
        if sequence_name == "Ref":
            sequence_name += "seq"
        st.markdown(f"""
            <style>
                .sequence-dashboard {{
                    font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                    margin-bottom: {nm_sizes['dashboard_margin']};
                }}
                .sequence-container {{
                    border: 1px solid #e1e5e9;
                    border-radius: {nm_sizes['container_border_radius']};
                    overflow: hidden;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.05);
                    background: white;
                    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                }}
                .sequence-container:hover {{
                    box-shadow: 0 4px 18px rgba(0,0,0,0.08);
                    transform: translateY(-1px);
                }}
                .sequence-header {{
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                    padding: {nm_sizes['header_padding']};
                    font-weight: 700;
                    font-size: {nm_sizes['header_font']};
                    display: flex;
                    gap: {nm_sizes['header_gap']};
                    align-items: center;
                }}
                .sequence-info-short {{
                    display: flex;
                    gap: 8px;
                    align-items: center;
                }}
                .sequence-length {{
                    font-size: {nm_sizes['length_font']};
                    font-weight: 900;
                    color: #fffae3;
                    background: rgba(255,255,255,0.22);
                    padding: {nm_sizes['length_padding']};
                    border-radius: {nm_sizes['length_border_radius']};
                    letter-spacing:{nm_sizes['length_letter_spacing']};
                    box-shadow: 0 0px 2px rgba(0,0,0,0.04);
                }}
                .no-motif-section {{
                    padding: {nm_sizes['section_padding']};
                    background: linear-gradient(135deg, #fff7f7 0%, #feebeb 100%);
                    border: 1px solid #feb2b2;
                    border-radius: {nm_sizes['section_border_radius']};
                    margin: {nm_sizes['section_margin']};
                    text-align: center;
                }}
                .no-motif-icon {{
                    font-size: {nm_sizes['icon_size']};
                    margin-bottom: 6px;
                }}
                .no-motif-msg {{
                    color: #c53030;
                    font-weight: 600;
                    font-size: {nm_sizes['message_font']};
                    margin-bottom: 6px;
                }}
                .no-motif-desc {{
                    color: #744210;
                    font-size: {nm_sizes['desc_font']};
                    opacity: {nm_sizes['desc_opacity']};
                }}
                .sequence-scroll-wrapper {{
                    width: 100%;
                    overflow-x: auto;
                    overflow-y: visible;
                    border-radius: 12px;
                    scrollbar-width: thin;
                    scrollbar-color: #cbd5e0 #f8fafc;
                }}
                .sequence-scroll-wrapper::-webkit-scrollbar {{
                    height: {scrollbar_sizes['height']};
                }}
                .sequence-scroll-wrapper::-webkit-scrollbar-thumb {{
                    background: #cbd5e0;
                }}
                .sequence-content-wrapper {{
                    display: inline-block;
                    overflow: visible;
                }}
                .sequence-content {{
                    padding: {nm_sizes['content_padding']};
                    white-space: nowrap;
                    overflow-x: visible;
                    overflow-y: hidden;
                    background: #f8fafc;
                    line-height: {nm_sizes['content_line_height']};
                    font-size: {nm_sizes['content_font']};
                    font-weight: 500;
                }}
                .sequence-scale {{
                    position: relative;
                    height: {nm_sizes['scale_height']};
                    padding: {nm_sizes['scale_padding']};
                    font-size: {nm_sizes['scale_font']};
                    color: #718096;
                    font-weight: 600;
                    background: #f8fafc;
                    white-space: nowrap;
                    overflow: visible;
                    min-height: {nm_sizes['scale_min_height']};
                }}
                .sequence-scale .scale-marker {{
                    position: absolute;
                    transform: translateX(-50%);
                    white-space: nowrap;
                    z-index: 10;
                }}
                .sequence-scale .scale-marker:first-child {{ 
                    transform: none; 
                    left: {nm_sizes['scale_marker_offset']} !important;
                }}
                .sequence-scale .scale-marker:last-child {{ 
                    transform: none; 
                    right: {nm_sizes['scale_marker_offset']} !important;
                }}
                .sequence-content::-webkit-scrollbar {{
                    height: {scrollbar_sizes['height']};
                }}
                .sequence-content::-webkit-scrollbar-thumb {{
                    background: #cbd5e0;
                    border-radius: 8px;
                }}
                .sequence-content::-webkit-scrollbar-track {{
                    background: #f8fafc;
                }}
            </style>
            
            <div class="sequence-dashboard">
                <div class="sequence-container">
                    <div class="sequence-header">
                        <div class="sequence-info-short">
                            <span>{sequence_name}</span>
                            <span class="sequence-length">{len(sequence)} base{'' if len(sequence)==1 else 's'}</span>
                        </div>
                    </div>
                    <div class="no-motif-section">
                        <div class="no-motif-icon">🔍</div>
                        <div class="no-motif-msg">No motifs detected in this region</div>
                        <div class="no-motif-desc">The sequence contains no recognizable motif patterns</div>
                    </div>
                    <div class="sequence-scroll-wrapper">
                        <div class="sequence-content-wrapper">
                            <div class="sequence-content" id="{sequence_content_id}">
                                <span style="color:#4a5568; letter-spacing:0.5px;">{sequence}</span>
                            </div>
                            <div class="sequence-scale" id="{sequence_scale_id}">
                                <span class="scale-marker" style="left:{nm_sizes['scale_marker_offset']};">0 bp</span>
                                <span class="scale-marker" style="left:25%;">{int(len(sequence) * 0.25)} bp</span>
                                <span class="scale-marker" style="left:50%;">{int(len(sequence) * 0.5)} bp</span>
                                <span class="scale-marker" style="left:75%;">{int(len(sequence) * 0.75)} bp</span>
                                <span class="scale-marker" style="right:{nm_sizes['scale_marker_offset']};">{len(sequence)} bp</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        <script>
            (function() {{
                const sequenceContentId = "{sequence_content_id}";
                const sequenceScaleId = "{sequence_scale_id}";

                function syncScaleWidth() {{
                    const sequenceContent = document.getElementById(sequenceContentId);
                    const sequenceScale = document.getElementById(sequenceScaleId);
                    const contentWrapper = sequenceContent ? sequenceContent.parentElement : null;

                    if (sequenceContent && sequenceScale && contentWrapper) {{
                        sequenceContent.offsetHeight; // force reflow

                        const contentWidth = sequenceContent.scrollWidth;
                        const scrollWrapper = contentWrapper.parentElement;

                        if (scrollWrapper) {{
                            const allWrappers = scrollWrapper.querySelectorAll('.sequence-content-wrapper');
                            let maxWidth = contentWidth;

                            allWrappers.forEach(wrapper => {{
                                const content = wrapper.querySelector('.sequence-content');
                                if (content) {{
                                    const width = content.scrollWidth;
                                    if (width > maxWidth) maxWidth = width;
                                }}
                            }});

                            allWrappers.forEach(wrapper => {{
                                const content = wrapper.querySelector('.sequence-content');
                                const scale = wrapper.querySelector('.sequence-scale');
                                if (content && scale) {{
                                    wrapper.style.width = maxWidth + 'px';
                                    scale.style.width = maxWidth + 'px';
                                    scale.style.minWidth = maxWidth + 'px';
                                }}
                            }});

                            const sequenceContainer = scrollWrapper.parentElement;
                            if (sequenceContainer && sequenceContainer.classList.contains('sequence-container')) {{
                                const allContainers = Array.from(document.querySelectorAll('.sequence-container'));
                                const sameLevelContainers = allContainers.filter(container => container.querySelector('.sequence-scroll-wrapper') !== null);

                                sameLevelContainers.forEach(container => {{
                                    const header = container.querySelector('.sequence-header');
                                    const statsBar = container.querySelector('.stats-bar');
                                    if (header) header.style.width = maxWidth + 'px';
                                    if (statsBar) statsBar.style.width = maxWidth + 'px';
                                }});
                            }} else {{
                                const parent = scrollWrapper.parentElement;
                                if (parent) {{
                                    const header = parent.querySelector('.sequence-header');
                                    const statsBar = parent.querySelector('.stats-bar');
                                    if (header) header.style.width = maxWidth + 'px';
                                    if (statsBar) statsBar.style.width = maxWidth + 'px';
                                }}
                            }}
                        }} else {{
                            contentWrapper.style.width = contentWidth + 'px';
                            sequenceScale.style.width = contentWidth + 'px';
                            sequenceScale.style.minWidth = contentWidth + 'px';

                            const sequenceContainer = sequenceContent.closest('.sequence-container');
                            if (sequenceContainer) {{
                                const header = sequenceContainer.querySelector('.sequence-header');
                                const statsBar = sequenceContainer.querySelector('.stats-bar');
                                if (header) header.style.width = contentWidth + 'px';
                                if (statsBar) statsBar.style.width = contentWidth + 'px';
                            }}
                        }}
                    }}
                }}

                document.addEventListener('DOMContentLoaded', syncScaleWidth);
                setTimeout(syncScaleWidth, 100);

                let resizeTimeout;
                window.addEventListener('resize', function() {{
                    clearTimeout(resizeTimeout);
                    resizeTimeout = setTimeout(syncScaleWidth, 50);
                }});

                if (typeof ResizeObserver !== 'undefined') {{
                    const resizeObserver = new ResizeObserver(syncScaleWidth);
                    const sequenceContent = document.getElementById(sequenceContentId);
                    if (sequenceContent) {{
                        setTimeout(() => resizeObserver.observe(sequenceContent), 200);
                    }}
                }}
            }})();
        </script>
        """, unsafe_allow_html=True)
        return

    # Process sequences with motifs
    ranges = parse_motif_range(spans)
    highlighted_sequence = ""
    previous_end = 0

    interruption_class = "interruption-segment motif-segment"
    interruption_style = "background: linear-gradient(135deg, #fc8181, #e53e3e); opacity: 0.85; border: 1px dashed rgba(255,255,255,0.5); color: #fff;"

    def calculate_gc_content(seq: str) -> float:
        if not seq:
            return 0.0
        gc_count = sum(1 for base in seq.upper() if base in {"G", "C"})
        return (gc_count / len(seq)) * 100

    total_motif_segments = len(ranges)

    # build the highlighted sequence
    for idx, (start, end) in enumerate(ranges):
        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        motif_name = motif_names[int(motif)]

        # add interruption if needed
        if start > previous_end:
            interruption_sequence = sequence[previous_end:start]
            int_start_pos_1based = previous_end + 1
            int_end_pos_1based = start
            interruption_title = html.escape(
                f"Interruption region | {int_start_pos_1based}-{int_end_pos_1based} | {len(interruption_sequence)} bp | {interruption_sequence}"
            )
            interruption_tooltip = (
                f"⚡ <strong>Interruption Region</strong>"
                f"<br/>📍 Position: {int_start_pos_1based}-{int_end_pos_1based}"
                f"<br/>📏 Length: {len(interruption_sequence)} bases"
                f"<br/>🧪 Sequence: <code>{interruption_sequence}</code>"
                f"<br/>🔍 Type: Non-motif sequence"
            )
            highlighted_sequence += (
                f"<span class='{interruption_class}' data-type='interruption' "
                f"style='{interruption_style}' "
                f"title=\"{interruption_title}\" "
                f"data-content=\"{interruption_tooltip}\">"
                f"<span class='interruption-text'>{interruption_sequence}</span>"
                f"</span>"
            )

        # add motif
        motif_sequence = sequence[start:end+1]
        motif_length = len(motif_sequence)

        motif_gc = calculate_gc_content(motif_sequence)

        # positions are 1-based for display
        start_pos_1based = start + 1
        end_pos_1based = end + 1
        
        motif_title = html.escape(
            f"{motif_name} | Position {start_pos_1based}-{end_pos_1based} | Length {motif_length} bp | Sequence {motif_sequence} | ID {motif} | Copy {idx + 1}/{total_motif_segments} | GC {motif_gc:.1f}%"
        )
        motif_tooltip = (
            f"🧬 <strong>{motif_name}</strong>"
            f"<br/>📍 Position: {start_pos_1based}-{end_pos_1based}"
            f"<br/>📏 Length: {motif_length} bases"
            f"<br/>🧪 Sequence: <code>{motif_sequence}</code>"
            f"<br/>🆔 Motif ID: {motif}"
            f"<br/>🔁 Copy index: {idx + 1} of {total_motif_segments}"
            f"<br/>🧮 GC content: {motif_gc:.1f}%"
        )
        highlighted_sequence += (
            f"<span class='motif-segment motif-{motif} motif-length-{motif_length}' data-motif='{motif}' "
            f"style='background-color:{color};' "
            f"title=\"{motif_title}\" "
            f"data-content=\"{motif_tooltip}\">"
            f"<span class='motif-text'>{motif_sequence}</span>"
            f"</span>"
        )
        previous_end = end + 1

    # add remaining sequence as interruption
    if previous_end < len(sequence):
        interruption_sequence = sequence[previous_end:]
        # final interruption positions
        final_int_start_pos_1based = previous_end + 1
        final_int_end_pos_1based = len(sequence)
        interruption_title = html.escape(
            f"Interruption region | {final_int_start_pos_1based}-{final_int_end_pos_1based} | {len(interruption_sequence)} bp | {interruption_sequence}"
        )
        interruption_tooltip = (
            f"⚡ <strong>Interruption Region</strong>"
            f"<br/>📍 Position: {final_int_start_pos_1based}-{final_int_end_pos_1based}"
            f"<br/>📏 Length: {len(interruption_sequence)} bases"
            f"<br/>🧪 Sequence: <code>{interruption_sequence}</code>"
            f"<br/>🔍 Type: Non-motif sequence"
        )
        highlighted_sequence += (
            f"<span class='{interruption_class}' data-type='interruption' "
            f"style='{interruption_style}' "
            f"title=\"{interruption_title}\" "
            f"data-content=\"{interruption_tooltip}\">"
            f"<span class='interruption-text'>{interruption_sequence}</span>"
            f"</span>"
        )

    if sequence_name == "Ref":
        sequence_name += "seq"

    # Calculate stats for motifs
    total_motifs = len(motif_ids)
    coverage = sum([end - start + 1 for start, end in ranges]) / len(sequence) * 100

    # Show supporting reads if available
    supporting_reads_html = ""
    if supporting_reads is not None:
        try:
            supporting_reads_val = int(supporting_reads)
            supporting_reads_html = (
                f"<div class='stat-item'>"
                f"    <span>Supporting reads:</span>"
                f"    <span class='stat-value'>{supporting_reads_val}</span>"
                f"</div>"
            )
        except Exception:
            supporting_reads_html = ""

    st.html(f"""
        <style>
            .sequence-dashboard {{
                font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                margin-bottom: {seq_sizes['dashboard_margin']};
            }}
            .sequence-container {{
                border: 1px solid #e1e5e9;
                border-radius: {seq_sizes['container_border_radius']};
                overflow: hidden;
                box-shadow: 0 2px 10px rgba(0,0,0,0.05);
                background: white;
            }}
            .sequence-header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: {seq_sizes['header_padding']};
                font-weight: 700;
                font-size: {seq_sizes['header_font']};
                display: flex;
                justify-content: space-between;
                align-items: center;
            }}
            .sequence-length {{
                font-size: {seq_sizes['length_font']};
                opacity: 0.92;
                font-weight: 500;
                background: rgba(255,255,255,0.15);
                padding: {seq_sizes['length_padding']};
                border-radius: {seq_sizes['length_border_radius']};
            }}
            .stats-bar {{
                display: flex;
                gap: {seq_sizes['stats_gap']};
                padding: {seq_sizes['stats_padding']};
                background: #f0f4ff;
                border-bottom: 1px solid #e1e8ff;
                font-size: {seq_sizes['stats_font']};
                color: #4a5568;
            }}
            .sequence-scroll-wrapper {{
                width: 100%;
                overflow-x: auto;
                overflow-y: visible;
                border-radius: 12px;
            }}
            .sequence-scroll-wrapper::-webkit-scrollbar {{
                height: {scrollbar_sizes['height']};
            }}
            .sequence-content-wrapper {{
                display: inline-block;
                overflow: visible;
            }}
            .sequence-content::-webkit-scrollbar {{
                height: {scrollbar_sizes['height']};
            }}
            .sequence-content {{
                padding: {seq_sizes['content_padding']};
                white-space: nowrap;
                overflow-x: visible;
                overflow-y: hidden;
                background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%);
                line-height: {seq_sizes['content_line_height']};
                font-size: {seq_sizes['content_font']};
                font-weight: 500;
                min-height: {seq_sizes['content_min_height']};
            }}
            .sequence-scale {{
                position: relative;
                height: {seq_sizes['scale_height']};
                padding: {seq_sizes['scale_padding']};
                font-size: {seq_sizes['scale_font']};
                color: #718096;
                font-weight: 600;
                background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%);
                white-space: nowrap;
                overflow: visible;
                min-height: {seq_sizes['scale_min_height']};
            }}
            .sequence-scale .scale-marker {{
                position: absolute;
                transform: translateX(-50%);
                white-space: nowrap;
                z-index: 10;
            }}
            .sequence-scale .scale-marker:first-child {{ 
                transform: none; 
                left: {seq_sizes['scale_marker_offset']} !important;
            }}
            .sequence-scale .scale-marker:last-child {{ 
                transform: none; 
                right: {seq_sizes['scale_marker_offset']} !important;
            }}
            .motif-segment {{
                display: inline-block;
                padding: {seq_sizes['motif_padding']};
                margin: 0;
                border-radius: {seq_sizes['motif_border_radius']};
                font-weight: 800;
                color: #ffffff;
                box-shadow: 0 4px 12px rgba(0,0,0,0.2);
                transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
                cursor: pointer;
                border: 3px solid rgba(255,255,255,0.5);
                position: relative;
                overflow: hidden;
                font-size: {seq_sizes['motif_font']};
                text-align: center;
                vertical-align: middle;
                text-shadow: 0 2px 4px rgba(0,0,0,0.4);
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
            }}
            .motif-segment:hover {{
                transform: translateY(-3px) scale(1.08);
                box-shadow: 0 12px 35px rgba(0,0,0,0.3);
                z-index: 15;
            }}
            .interruption-segment {{
                color: #fff !important;
                background: linear-gradient(135deg, #fc8181, #e53e3e) !important;
                opacity: 0.95 !important;
                border: 3px dashed rgba(255,255,255,0.7) !important;
                border-radius: {seq_sizes['interruption_border_radius']} !important;
                margin: 0 !important;
                padding: {seq_sizes['interruption_padding']} !important;
                box-shadow: 0 4px 12px rgba(229, 62, 62, 0.3) !important;
            }}
            .interruption-segment:hover {{
                opacity: 1 !important;
                transform: translateY(-3px) scale(1.08) !important;
                box-shadow: 0 12px 35px rgba(229, 62, 62, 0.4) !important;
            }}
            .tooltip {{
                position: fixed;
                background: linear-gradient(135deg, rgba(15, 23, 42, 0.98), rgba(30, 41, 59, 0.98));
                color: white;
                padding: {tooltip_sizes['padding']};
                border-radius: {tooltip_sizes['border_radius']};
                font-size: {tooltip_sizes['font']};
                pointer-events: none;
                opacity: 0;
                transform: translateY(15px) scale(0.95);
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                z-index: 1000;
                backdrop-filter: blur(20px);
                border: 2px solid rgba(255,255,255,0.2);
                max-width: {tooltip_sizes['max_width']};
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
                box-shadow: 0 10px 40px rgba(0,0,0,0.3);
                line-height: 1.6;
            }}
            .tooltip.show {{
                opacity: 1;
                transform: translateY(0) scale(1);
            }}
            .tooltip strong {{
                color: #60a5fa;
            }}
            .tooltip code {{
                background: rgba(255,255,255,0.1);
                padding: {tooltip_sizes['code_padding']};
                border-radius: 4px;
                font-size: {tooltip_sizes['code_font']};
                color: #fbbf24;
            }}
        </style>

        <div class="sequence-dashboard">
            <div class="sequence-container">
                <div class="sequence-header">
                    <span>{sequence_name}</span>
                    <span class="sequence-length">{len(sequence)} bases</span>
                </div>
                <div class="stats-bar">
                    <div class="stat-item">
                        <span>Copy number:</span>
                        <span class="stat-value">{total_motifs}</span>
                    </div>
                    <div class="stat-item">
                        <span>Coverage:</span>
                        <span class="stat-value">{coverage:.1f}%</span>
                    </div>
                    {supporting_reads_html}
                </div>
                <div class="sequence-scroll-wrapper">
                    <div class="sequence-content-wrapper">
                        <div class="sequence-content" id="{sequence_content_id}">
                            {highlighted_sequence}
                        </div>
                        <div class="sequence-scale" id="{sequence_scale_id}">
                            <span class="scale-marker" style="left:{seq_sizes['scale_marker_offset']};">0 bp</span>
                            <span class="scale-marker" style="left:25%;">{int(len(sequence) * 0.25)} bp</span>
                            <span class="scale-marker" style="left:50%;">{int(len(sequence) * 0.5)} bp</span>
                            <span class="scale-marker" style="left:75%;">{int(len(sequence) * 0.75)} bp</span>
                            <span class="scale-marker" style="right:{seq_sizes['scale_marker_offset']};">{len(sequence)} bp</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        <script>
            (function() {{
                const sequenceContentId = "{sequence_content_id}";
                const sequenceScaleId = "{sequence_scale_id}";

                var currentTooltip = window.__sequenceTooltipEl || null;

                function createTooltip() {{
                    if (!currentTooltip) {{
                        currentTooltip = document.createElement('div');
                        currentTooltip.className = 'tooltip';
                        document.body.appendChild(currentTooltip);
                        window.__sequenceTooltipEl = currentTooltip;
                    }}
                    return currentTooltip;
                }}

                function showTooltip(content, x, y) {{
                    const tooltip = createTooltip();
                    tooltip.innerHTML = content;
                    tooltip.style.left = x + 'px';
                    tooltip.style.top = y + 'px';
                    tooltip.classList.add('show');
                }}

                function hideTooltip() {{
                    if (currentTooltip) {{
                        currentTooltip.classList.remove('show');
                    }}
                }}

                function initSequenceTooltips() {{
                    const container = document.getElementById(sequenceContentId);
                    if (!container) return;

                    const segments = container.querySelectorAll('[data-content]');

                    segments.forEach(segment => {{
                        segment.addEventListener('mouseenter', function(e) {{
                            const content = this.getAttribute('data-content');
                            showTooltip(content, e.pageX + 15, e.pageY + 15);
                        }});

                        segment.addEventListener('mouseleave', hideTooltip);

                        segment.addEventListener('mousemove', function(e) {{
                            if (currentTooltip && currentTooltip.classList.contains('show')) {{
                                currentTooltip.style.left = (e.pageX + 15) + 'px';
                                currentTooltip.style.top = (e.pageY + 15) + 'px';
                            }}
                        }});
                    }});
                }}

                function syncScaleWidth() {{
                    const sequenceContent = document.getElementById(sequenceContentId);
                    const sequenceScale = document.getElementById(sequenceScaleId);
                    const contentWrapper = sequenceContent ? sequenceContent.parentElement : null;

                    if (sequenceContent && sequenceScale && contentWrapper) {{
                        sequenceContent.offsetHeight;

                        const contentWidth = sequenceContent.scrollWidth;
                        const scrollWrapper = contentWrapper.parentElement;

                        if (scrollWrapper) {{
                            const allWrappers = scrollWrapper.querySelectorAll('.sequence-content-wrapper');
                            let maxWidth = contentWidth;

                            allWrappers.forEach(wrapper => {{
                                const content = wrapper.querySelector('.sequence-content');
                                if (content) {{
                                    const width = content.scrollWidth;
                                    if (width > maxWidth) maxWidth = width;
                                }}
                            }});

                            allWrappers.forEach(wrapper => {{
                                const content = wrapper.querySelector('.sequence-content');
                                const scale = wrapper.querySelector('.sequence-scale');
                                if (content && scale) {{
                                    wrapper.style.width = maxWidth + 'px';
                                    scale.style.width = maxWidth + 'px';
                                    scale.style.minWidth = maxWidth + 'px';
                                }}
                            }});

                            const sequenceContainer = scrollWrapper.parentElement;
                            if (sequenceContainer && sequenceContainer.classList.contains('sequence-container')) {{
                                const allContainers = Array.from(document.querySelectorAll('.sequence-container'));
                                const sameLevelContainers = allContainers.filter(container => container.querySelector('.sequence-scroll-wrapper') !== null);

                                sameLevelContainers.forEach(container => {{
                                    const header = container.querySelector('.sequence-header');
                                    const statsBar = container.querySelector('.stats-bar');
                                    if (header) header.style.width = maxWidth + 'px';
                                    if (statsBar) statsBar.style.width = maxWidth + 'px';
                                }});
                            }} else {{
                                const parent = scrollWrapper.parentElement;
                                if (parent) {{
                                    const header = parent.querySelector('.sequence-header');
                                    const statsBar = parent.querySelector('.stats-bar');
                                    if (header) header.style.width = maxWidth + 'px';
                                    if (statsBar) statsBar.style.width = maxWidth + 'px';
                                }}
                            }}
                        }} else {{
                            contentWrapper.style.width = contentWidth + 'px';
                            sequenceScale.style.width = contentWidth + 'px';
                            sequenceScale.style.minWidth = contentWidth + 'px';

                            const sequenceContainer = sequenceContent.closest('.sequence-container');
                            if (sequenceContainer) {{
                                const header = sequenceContainer.querySelector('.sequence-header');
                                const statsBar = sequenceContainer.querySelector('.stats-bar');
                                if (header) header.style.width = contentWidth + 'px';
                                if (statsBar) statsBar.style.width = contentWidth + 'px';
                            }}
                        }}
                    }}
                }}

                document.addEventListener('DOMContentLoaded', function() {{
                    initSequenceTooltips();
                    syncScaleWidth();
                }});

                setTimeout(function() {{
                    initSequenceTooltips();
                    syncScaleWidth();
                }}, 100);

                let resizeTimeout;
                window.addEventListener('resize', function() {{
                    clearTimeout(resizeTimeout);
                    resizeTimeout = setTimeout(syncScaleWidth, 50);
                }});

                if (typeof ResizeObserver !== 'undefined') {{
                    const resizeObserver = new ResizeObserver(function() {{
                        syncScaleWidth();
                    }});
                    const sequenceContent = document.getElementById(sequenceContentId);
                    if (sequenceContent) {{
                        setTimeout(() => resizeObserver.observe(sequenceContent), 200);
                    }}
                }}
            }})();
        </script>
    """)



def display_motifs_as_bars(sequence_name, motif_colors, motif_ids, spans, sequence, motif_names, supporting_reads=None):
    """
    This function draws motif bars for each motif in the sequence. 
    Motif information (name, count, coverage, sequence) is shown in a stylish tooltip on hover, not on the bar itself.
    This revision uses classic mouseover/mouseout events for tooltips and makes sure tooltips are attached per motif bar/interruption element.
    """

    sequence_length = len(sequence)
    ranges = parse_motif_range(spans)

    if not isinstance(motif_names, list):
        motif_names = [motif_names]

    # Count all motif occurrences (by id) in the sequence
    motif_counter = Counter(motif_ids)
    # Compute covered bases by motif
    motif_coverage = {}
    for idx, motif in enumerate(motif_ids):
        span = ranges[idx]
        motif_coverage[motif] = motif_coverage.get(motif, 0) + (span[1] - span[0] + 1)
    # Precompute tooltip info for each motif id
    motif_tooltip_map = {}
    for motif in set(motif_ids):
        # Name, count, total bases for this motif
        name = motif_names[int(motif)]
        count = motif_counter[motif]
        coverage_bases = motif_coverage[motif]
        coverage_percent = coverage_bases / sequence_length * 100
        motif_tooltip_map[motif] = (
            f"<b>{name}</b><br>"
            f"Count: <b>{count}</b><br>"
            f"Coverage: <b>{coverage_bases} bp</b> ({coverage_percent:.2f}%)"
        )

    motif_bar_htmls = []
    previous_end = 0
    gap = 0.3   # minimum gap between bars
    total_motifs = len(motif_ids)
    overall_coverage = sum([end - start + 1 for start, end in ranges]) / len(sequence) * 100

    for idx, (start, end) in enumerate(ranges):
        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        span_length = end - start + 1

        if start >= 0 and end <= sequence_length:
            if start > previous_end:
                interruption_width = (start - previous_end) / sequence_length * 100
                interruption_start = previous_end / sequence_length * 100
                motif_bar_htmls.append(
                    f"<div class='interruption-bar' data-motif='interruption' "
                    f"data-tooltip='Interruption: {sequence[previous_end:start]} ({start-previous_end} bases)' "
                    f"style='left:{interruption_start}%; width:{interruption_width}%;'>"
                    f"<div class='bar-pattern'></div>"
                    f"</div>"
                )

            relative_width = max((span_length / sequence_length) * 100 - gap, 0.5)
            relative_start = (start / sequence_length) * 100

            display_content = ""
            tooltip_html = (
                motif_tooltip_map[motif]
                + f"<br>Sequence region: <span style='font-family:monospace;'>{sequence[start:end+1]}</span><br>Start: <b>{start+1}</b> End: <b>{end+1}</b> ({span_length} bases)"
            )

            motif_bar_htmls.append(
                f"<div class='motif-bar motif-{motif}' data-motif='{motif}' "
                f"data-tooltip=\"{tooltip_html}\" "
                f"style='left:{relative_start}%; width:{relative_width}%; background-color:{color};'>"
                f"<div class='bar-glow'></div>"
                f"{display_content}"
                f"</div>"
            )

            previous_end = end + 1

    if previous_end < sequence_length:
        interruption_width = (sequence_length - previous_end) / sequence_length * 100
        interruption_start = previous_end / sequence_length * 100
        motif_bar_htmls.append(
            f"<div class='interruption-bar' data-motif='interruption' "
            f"data-tooltip='Interruption: {sequence[previous_end:]} ({sequence_length-previous_end} bases)' "
            f"style='left:{interruption_start}%; width:{interruption_width}%;'>"
            f"<div class='bar-pattern'></div>"
            f"</div>"
        )

    # Prepare supporting reads count html
    supporting_reads_html = ""
    if supporting_reads is not None:
        try:
            supporting_reads_val = int(supporting_reads)
            supporting_reads_html = (
                f"<div class='stat-item'>"
                f"    <span>Supporting reads:</span>"
                f"    <span class='stat-value'>{supporting_reads_val}</span>"
                f"</div>"
            )
        except Exception:
            supporting_reads_html = ""

    bar_container = f"""
    <div class="bar-visualization">
        <div class="sequence-bar-container" id="bar-container-{sequence_name}">
            <div class="bar-track">
                {''.join(motif_bar_htmls)}
            </div>
            <div class="bar-scale">
                <span>0</span>
                <span style="left:25%">25%</span>
                <span style="left:50%">50%</span>
                <span style="left:75%">75%</span>
                <span style="right:0">100%</span>
            </div>
        </div>
        <div id="tooltip-{sequence_name}" class="tooltip"></div>
    </div>
    """

    if sequence_name == "Ref":
        sequence_name += "seq"

    st.html(f"""
        <style>
            .sequence-dashboard {{
                font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                margin-bottom: 18px;
            }}
            .bar-visualization {{
                font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                margin-bottom: 18px;
            }}
            .sequence-container {{
                border: 1px solid #e1e5e9;
                border-radius: 12px;
                overflow: hidden;
                box-shadow: 0 2px 10px rgba(0,0,0,0.05);
                background: white;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            }}
            .sequence-container:hover {{
                box-shadow: 0 4px 18px rgba(0,0,0,0.08);
                transform: translateY(-1px);
            }}
            .sequence-header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 9px 12px;
                font-weight: 700;
                font-size: 14px;
                display: flex;
                justify-content: space-between;
                align-items: center;
            }}
            .sequence-length {{
                font-size: 14px;
                opacity: 0.92;
                font-weight: 500;
                background: rgba(255,255,255,0.15);
                padding: 2px 9px;
                border-radius: 14px;
            }}
            .stats-bar {{
                display: flex;
                gap: 11px;
                padding: 8px 12px;
                background: #f0f4ff;
                border-bottom: 1px solid #e1e8ff;
                font-size: 14px;
                color: #4a5568;
            }}
            .stat-item {{
                display: flex;
                align-items: center;
                gap: 4px;
                font-weight: 500;
            }}
            .stat-value {{
                background: white;
                padding: 1px 6px;
                border-radius: 8px;
                border: 1px solid #cbd5e0;
                font-weight: 600;
                color: #2d3748;
                font-size: 13px;
            }}
            .sequence-bar-container {{
                position: relative;
                height: 85px;
                background: #f8fafc;
                border-radius: 8px;
                padding: 10px 10px 15px 10px;
                margin: 0;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            }}
            .bar-track {{
                position: relative;
                height: 40px;
                background: #edf2f7;
                border-radius: 12px;
                overflow: hidden;
                border: 1px solid #cbd5e0;
            }}
            .motif-bar {{
                position: absolute;
                height: 36px;
                top: 2px;
                border-radius: 10px;
                border: 2px solid white;
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                cursor: pointer;
                display: flex;
                align-items: center;
                justify-content: center;
                font-weight: 700;
                font-size: 11px;
                color: white;
                text-shadow: 0 1px 2px rgba(0,0,0,0.3);
                overflow: hidden;
                min-width: 4px;
            }}
            .motif-bar:hover {{
                transform: scale(1.05) translateY(-2px);
                box-shadow: 0 8px 25px rgba(0,0,0,0.25);
                z-index: 10;
            }}
            .bar-glow {{
                position: absolute;
                top: 0;
                left: -100%;
                width: 100%;
                height: 100%;
                background: linear-gradient(90deg, transparent, rgba(255,255,255,0.4), transparent);
                transition: left 0.6s ease;
            }}
            .motif-bar:hover .bar-glow {{
                left: 100%;
            }}
            .bar-label {{
                font-size: 11px;
                font-weight: 700;
                color: white;
                text-shadow: 0 1px 2px rgba(0,0,0,0.5);
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
                max-width: 100%;
                padding: 0 2px;
                text-align: center;
            }}
            .bar-dot {{
                font-size: 16px;
                color: white;
                text-shadow: 0 1px 2px rgba(0,0,0,0.5);
            }}
            .interruption-bar {{
                position: absolute;
                height: 36px;
                top: 2px;
                background: #fed7d7;
                border: 2px dashed #e53e3e;
                border-radius: 8px;
                cursor: pointer;
                transition: all 0.3s ease;
                overflow: hidden;
            }}
            .interruption-bar:hover {{
                background: #feebeb;
                transform: scaleY(1.1);
            }}
            .bar-pattern {{
                width: 100%;
                height: 100%;
                background: repeating-linear-gradient(
                    45deg,
                    transparent,
                    transparent 5px,
                    rgba(229, 62, 62, 0.1) 5px,
                    rgba(229, 62, 62, 0.1) 10px
                );
            }}
            .bar-scale {{
                position: relative;
                height: 30px;
                margin-top: 12px;
                padding-top: 5px;
                font-size: 14px;
                color: #718096;
                font-weight: 600;
            }}
            .bar-scale span {{
                position: absolute;
                transform: translateX(-50%);
            }}
            .bar-scale span:first-child {{ left: 0; transform: none; }}
            .bar-scale span:last-child {{ left: auto; right: 0; transform: none; }}
            .tooltip {{
                position: fixed;
                background: rgba(45, 55, 72, 0.97);
                color: white;
                padding: 12px 16px;
                border-radius: 10px;
                font-size: 13px;
                pointer-events: none;
                opacity: 0;
                transform: translateY(10px);
                transition: all 0.2s ease;
                z-index: 1000;
                backdrop-filter: blur(10px);
                border: 1px solid rgba(255,255,255,0.07);
                max-width: 400px;
                font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                box-shadow: 0 6px 30px rgba(0,0,0,0.18);
                letter-spacing: 0.01em;
                line-height: 1.6;
            }}
            .tooltip.show {{
                opacity: 1;
                transform: translateY(0) scale(1);
            }}
            .highlighted {{
                transform: scale(1.1);
                box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.5);
                z-index: 5;
            }}
        </style>

        <div class="sequence-dashboard">
            <div class="sequence-container">
                <div class="sequence-header">
                    <span>{sequence_name}</span>
                    <span class="sequence-length">{sequence_length} bases</span>
                </div>
                <div class="stats-bar">
                    <div class="stat-item">
                        <span>Copy number:</span>
                        <span class="stat-value">{total_motifs}</span>
                    </div>
                    <div class="stat-item">
                        <span>Coverage:</span>
                        <span class="stat-value">{overall_coverage:.1f}%</span>
                    </div>
                    {supporting_reads_html}
                </div>
                {bar_container}
            </div>
        </div>

        <script>
            function initializeBarVisualization(sequenceName) {{
                // Wait for DOM to be ready
                setTimeout(() => {{
                    const container = document.getElementById('bar-container-' + sequenceName);
                    const tooltip = document.getElementById('tooltip-' + sequenceName);
                    
                    if (!container || !tooltip) {{
                        console.log('Container or tooltip not found for:', sequenceName);
                        return;
                    }}

                    // Helper to get position
                    function getOffset(evt) {{
                        if ('touches' in evt && evt.touches.length > 0) {{
                            return {{ x: evt.touches[0].clientX, y: evt.touches[0].clientY }};
                        }} else {{
                            return {{ x: evt.clientX, y: evt.clientY }};
                        }}
                    }}

                    // Attach events to motif/interruption bars
                    const bars = container.querySelectorAll('.motif-bar, .interruption-bar');
                    console.log('Found bars:', bars.length);
                    
                    bars.forEach(bar => {{
                        const tooltipContent = bar.getAttribute('data-tooltip');
                        console.log('Bar tooltip content:', tooltipContent);
                        
                        bar.addEventListener('mouseenter', function(event) {{
                            if (tooltipContent) {{
                                tooltip.innerHTML = tooltipContent;
                                let coords = getOffset(event);
                                tooltip.style.left = (coords.x + 15) + 'px';
                                tooltip.style.top = (coords.y + 15) + 'px';
                                tooltip.classList.add('show');
                                console.log('Tooltip shown:', tooltipContent);
                            }}
                        }});
                        
                        bar.addEventListener('mousemove', function(event) {{
                            if (tooltip.classList.contains('show')) {{
                                let coords = getOffset(event);
                                tooltip.style.left = (coords.x + 15) + 'px';
                                tooltip.style.top = (coords.y + 15) + 'px';
                            }}
                        }});
                        
                        bar.addEventListener('mouseleave', function() {{
                            tooltip.classList.remove('show');
                        }});
                    }});
                }}, 100);
            }}

            initializeBarVisualization('{sequence_name}');
        </script>
    """)



def plot_motif_bar(motif_count, motif_names, motif_colors=None, sequence_name=""):
    motif_labels = []
    motif_counts = []
    motif_ids = []
    
    for label, value in sorted(motif_count.items()):
        motif_name = motif_names[int(label)]
        if motif_name:
            motif_labels.append(motif_name)
            motif_counts.append(value)
            motif_ids.append(int(label))
    
    data = {
        'Motif': motif_labels,
        'Count': motif_counts,
        'Motif_ID': motif_ids
    }
    df = pd.DataFrame(data)
    
    color_list = [motif_colors[motif_id] for motif_id in motif_ids] if motif_colors else None

    # Create interactive bar chart
    bar_chart = alt.Chart(df).mark_bar(
        cornerRadius=8,
        stroke='white',
        strokeWidth=2
    ).encode(
        y=alt.Y('Count:Q', 
            title='Occurrences',
            axis=alt.Axis(
                labelColor='#4B5563',
                titleColor='#374151',
                labelFontWeight='bold',
                titleFontWeight='bold',
                tickMinStep=1,              
                format='d'                 
            )),
        x=alt.X('Motif:N', 
            sort='-y',
            title='',
            axis=alt.Axis(
                labelColor='#4B5563',
                labelFontWeight='bold'
            )),
        color=alt.Color('Motif:N',
                    scale=alt.Scale(domain=motif_labels, range=color_list),
                    legend=None),
        tooltip=['Motif', 'Count']
    ).properties(
        width=400,
        height=300,
        title=alt.TitleParams(
            text=f'🧬 Motif Occurrences - {sequence_name}',
            fontSize=16,
            fontWeight='bold',
            color='#1F2937'
        )
    ).configure_view(
        strokeWidth=0,
        fill='rgba(255,255,255,0.9)'
    )

    # Plot the bar chart
    bar_chart = bar_chart.configure_axisX(labelAngle=0).configure_legend(disable=True)
    st.altair_chart(bar_chart, use_container_width=True)

def interpret_genotype(gt):
    """
    Interpret genotype string and return meaningful description.
    
    Args:
        gt (str): Genotype string like "0/0", "0/1", "1/2", "0", "1", etc.
    
    Returns:
        dict: Interpretation with description, colors, and icons
    """
    if not gt or gt in ["./.", ".", ""]:
        return {
            'description': 'No genotype called',
            'interpretation': 'Missing data',
            'color': '#9CA3AF',
            'bg_color': '#F3F4F6',
            'icon': '❓',
            'status': 'unknown'
        }
    
    # Clean and normalize the input
    gt_str = str(gt).strip()
    # For assembly mode, genotype may be a single allele ("0", "1") or a tuple/list
    if st.session_state.get("cohort_mode", "") == "assembly":
        if isinstance(gt, (list, tuple)):
            alleles = [str(a) for a in gt]
        else:
            # Handle "0" or "1", or "0/1" 
            if "/" in gt_str:
                alleles = gt_str.split("/")
            elif "|" in gt_str:
                alleles = gt_str.split("|")
            else:
                alleles = [gt_str]
    else:
        
        if "/" in gt_str:
            alleles = gt_str.split("/")
        elif "|" in gt_str:
            alleles = gt_str.split("|")
        else:
            alleles = [gt_str]

    # Remove empty alleles 
    alleles = [a for a in alleles if a not in [".", ""]]

    # Accept single-haplotype genotypes as well 
    if len(alleles) == 1:
        allele = alleles[0]
        if allele == "0":
            return {
                'description': 'Hemizygous Reference ',
                'interpretation': 'Only one haplotype; matches reference',
                'color': '#10B981',
                'bg_color': '#D1FAE5',
                'icon': '🟢',
                'status': 'hemizygous_ref'
            }
        else:
            return {
                'description': f'Hemizygous Alternative',
                'interpretation': f'Only one haplotype; differs from reference (allele {allele})',
                'color': '#F59E0B',
                'bg_color': '#FEF3C7',
                'icon': '🟡',
                'status': 'hemizygous_alt'
            }
    elif len(alleles) != 2:
        return {
            'description': f'Invalid genotype: {gt}',
            'interpretation': 'Malformed',
            'color': '#EF4444',
            'bg_color': '#FEE2E2',
            'icon': '⚠️',
            'status': 'error'
        }

    allele1, allele2 = alleles
    
    # Handle different genotype patterns
    if allele1 == allele2:
        if allele1 == "0":
            return {
                'description': 'Homozygous Reference ',
                'interpretation': 'Both haplotypes identical to reference',
                'color': '#10B981',
                'bg_color': '#D1FAE5',
                'icon': '🟢',
                'status': 'homozygous_ref'
            }
        else:
            return {
                'description': f'Homozygous Alternative ',
                'interpretation': f'Both haplotypes different from reference (allele {allele1})',
                'color': '#F59E0B',
                'bg_color': '#FEF3C7',
                'icon': '🟡',
                'status': 'homozygous_alt'
            }
    else:
        if allele1 == "0" or allele2 == "0":
            return {
                'description': f'Heterozygous Reference/Alternative ',
                'interpretation': 'One haplotype like reference, one different',
                'color': '#3B82F6',
                'bg_color': '#DBEAFE',
                'icon': '🔵',
                'status': 'heterozygous_ref_alt'
            }
        else:
            return {
                'description': f'Heterozygous Alternative ',
                'interpretation': f'Both haplotypes different from reference (alleles {allele1} and {allele2})',
                'color': '#8B5CF6',
                'bg_color': '#EDE9FE',
                'icon': '🟣',
                'status': 'heterozygous_alt_alt'
            }


def display_genotype_card(gt, sample_name="Sample", show_details=True):
    """
    Display genotype information in a compact card
    
    Args:
        gt (str): Genotype string
        sample_name (str): Name of the sample
        show_details (bool): Whether to show detailed interpretation
    """
    interpretation = interpret_genotype(gt)

    st.html(f"""
        <style>
            .genotype-card {{
                background: linear-gradient(135deg, {interpretation['bg_color']} 0%, rgba(255,255,255,0.92) 100%);
                border: 1px solid {interpretation['color']};
                border-radius: 4px;
                padding: 0.5px 3px; /* Minimally reduce padding */
                margin: 0.5px 0;
                box-shadow: 0 0.5px 1.5px rgba(0,0,0,0.03);
                font-size: 14px !important;
                transition: all 0.12s;
                position: relative;
                overflow: hidden;
                min-height: 0 !important;
                line-height: 1.03;
            }}
            .genotype-card:hover {{
                transform: translateY(-0.5px) scale(1.008);
                box-shadow: 0 1px 4px rgba(0,0,0,0.065);
            }}
            .genotype-content {{
                display: flex;
                align-items: center;
                gap: 1px; 
                flex-wrap: wrap;
                font-size: 14px !important;
                min-height: 0;
                line-height: 1.07;
            }}
            .genotype-icon {{
                font-size: 14px;
            }}
            .genotype-title {{
                font-size: 14px !important;
                font-weight: 600;
                color: {interpretation['color']};
                margin: 0 0.5px 0 0;
            }}
            .genotype-value {{
                background: {interpretation['color']};
                color: white;
                padding: 0px 3px;
                border-radius: 4px;
                font-size: 14px !important;
                font-weight: 700;
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
                letter-spacing: 0.07px;
                height: 16px;
                display: inline-flex;
                align-items: center;
            }}
            .genotype-stats {{
                display: flex;
                gap: 1.5px;
                margin-left: 7px;
                    font-size: 14px !important;
            }}
            .stat-item {{
                font-size: 7px !important;
                color: #6B7280;
                font-weight: 500;
                line-height: 1;
            }}
            .stat-value {{
                color: {interpretation['color']};
                font-weight: 600;
                font-size: 14px !important;
                line-height: 1;
            }}
        </style>
        
        <div class="genotype-card">
            <div class="genotype-content">
                <span class="genotype-icon">{interpretation['icon']}</span>
                <span class="genotype-title">{sample_name}:</span>
                <span class="genotype-value">{gt}</span>
                <div class="genotype-stats">
                    <span class="stat-item"><span class="stat-value">{interpretation['status'].replace('_', ' ').title()}</span></span>
                    <span class="stat-item">•</span>
                    <span class="stat-item"><span class="stat-value">{'Diploid' if '/' in gt else 'Unknown'}</span></span>
                </div>
            </div>
        </div>
    """)


def display_genotype_badge(gt, size="medium"):
    """
    Display genotype as a compact badge
    
    Args:
        gt (str): Genotype string
        size (str): Size of badge ("small", "medium", "large")
    """
    interpretation = interpret_genotype(gt)
    
    size_classes = {
        "small": "padding: 5px 10px; font-size: 14px; border-radius: 9px;",
        "medium": "padding: 7px 14px; font-size: 16px; border-radius: 11px;",
        "large": "padding: 9px 18px; font-size: 18px; border-radius: 13px;"
    }
    
    st.html(f"""
        <style>
            .genotype-badge {{
                display: inline-block;
                background: {interpretation['bg_color']};
                color: {interpretation['color']};
                border: 2px solid {interpretation['color']};
                font-weight: 700;
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
                letter-spacing: 0.5px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                transition: all 0.2s ease;
                cursor: pointer;
                {size_classes[size]}
            }}
            .genotype-badge:hover {{
                transform: translateY(-1px);
                box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            }}
        </style>
        
        <span class="genotype-badge" title="{interpretation['description']}: {interpretation['interpretation']}">
            {interpretation['icon']} {gt}
        </span>
    """)


def _compute_genotype_cache(genotypes_dict):
    """
    Compute and cache all genotype-related data once.
    Returns a dictionary with all pre-computed data.
    """
    samples = list(genotypes_dict.keys())
    unique_genotypes = list(set(genotypes_dict.values()))

    # Color mapping and genotype data
    genotype_colors = {}
    genotype_data = {}
    for gt in unique_genotypes:
        interpretation = interpret_genotype(gt)
        genotype_colors[gt] = interpretation['color']
        genotype_data[gt] = {
            'color': interpretation['color'],
            'bg_color': interpretation['bg_color'],
            'icon': interpretation['icon'],
            'status': interpretation['status'],
            'description': interpretation['description']
        }

    # Group samples by genotype
    genotype_groups = {}
    for sample, gt in genotypes_dict.items():
        if gt not in genotype_groups:
            genotype_groups[gt] = []
        genotype_groups[gt].append(sample)
    
    # Sort groups by count (descending)
    sorted_groups = sorted(genotype_groups.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Pre-compute sorted genotype list for filter dropdown (sorted by count descending)
    sorted_genotypes = sorted(unique_genotypes, key=lambda x: len(genotype_groups.get(x, [])), reverse=True)
    filter_options = ["All"] + sorted_genotypes
    
    # Pre-compute summary statistics
    homozygous_count = 0
    heterozygous_count = 0
    missing_count = 0
    for gt in genotypes_dict.values():
        gt_info = genotype_data.get(gt, {})
        status = gt_info.get('status', '')
        if status.startswith('homozygous'):
            homozygous_count += 1
        elif status.startswith('heterozygous'):
            heterozygous_count += 1
        elif status == 'unknown':
            missing_count += 1
    
    return {
        'samples': samples,
        'unique_genotypes': unique_genotypes,
        'filter_options': filter_options,  # Pre-computed sorted list with "All"
        'genotype_colors': genotype_colors,
        'genotype_data': genotype_data,
        'genotype_groups': genotype_groups,
        'sorted_groups': sorted_groups,
        'homozygous_count': homozygous_count,
        'heterozygous_count': heterozygous_count,
        'missing_count': missing_count
    }


def create_genotype_comparison_matrix(genotypes_dict):
    """
    Create a visual comparison matrix of genotypes across samples with multiple view modes
    Args:
        genotypes_dict (dict): Dictionary with sample names as keys and genotypes as values
    """

    st.markdown("###  Genotype Comparison Matrix")

    # Initialize or update cache in session state
    cache_key = 'genotype_matrix_cache'
    cached_dict_key = 'genotype_matrix_cached_dict'
    
    # Check if we need to recompute by comparing dictionaries directly
    if (cache_key not in st.session_state or 
        cached_dict_key not in st.session_state or 
        st.session_state[cached_dict_key] != genotypes_dict):
        # Compute and cache all data
        st.session_state[cache_key] = _compute_genotype_cache(genotypes_dict)
        st.session_state[cached_dict_key] = genotypes_dict.copy()
    
    # Get cached data
    cache = st.session_state[cache_key]
    samples = cache['samples']
    unique_genotypes = cache['unique_genotypes']
    genotype_data = cache['genotype_data']
    genotype_groups = cache['genotype_groups']
    sorted_groups = cache['sorted_groups']

    # Use tabs instead of radio - tabs only render active content, much faster!
    tab1, tab2 = st.tabs(["👨‍👩‍👧‍👦 Grouped View", "🗂️ Grid Cards"])
    
    # Use all samples - no filtering
    filtered_dict = genotypes_dict
    
    with tab1:
        _display_grouped_view(filtered_dict, sorted_groups, genotype_data, genotype_groups)
    
    with tab2:
        _display_grid_cards(filtered_dict, genotype_data)

    # Add summary statistics - use cached values
    st.markdown("#### 📊 Genotype Summary")
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.metric("Total Samples", len(samples))
    with col2:
        st.metric("Unique Genotypes", len(unique_genotypes))
    with col3:
        st.metric("Homozygous", cache['homozygous_count'])
    with col4:
        st.metric("Heterozygous", cache['heterozygous_count'])
    with col5:
        st.metric("Missing", cache['missing_count'])


def _display_grouped_view(filtered_dict, sorted_groups, genotype_data, genotype_groups):
    """Display genotypes grouped by genotype type with collapsible sections"""
    
    # Since we're not filtering, filtered_dict == genotypes_dict, so use sorted_groups directly
    for gt, group_samples in sorted_groups:
        if not group_samples:
            continue
        
        gt_info = genotype_data[gt]  # Use cached data
        
        # Create expander for each genotype group
        with st.expander(
            f"{gt_info['icon']} **{gt}** - {gt_info['description']} ({len(group_samples)} samples)",
            expanded=len(sorted_groups) <= 3  # Auto-expand if few groups
        ):
            # Display all samples in compact grid - no pagination
            cols = st.columns(min(6, len(group_samples)))
            for idx, sample in enumerate(group_samples):
                with cols[idx % len(cols)]:
                    st.markdown(f"""
                        <div style="
                            background: {gt_info['bg_color']};
                            border: 2px solid {gt_info['color']};
                            border-radius: 8px;
                            padding: 8px;
                            text-align: center;
                            margin-bottom: 8px;
                        ">
                            <div style="font-size: 14px; font-weight: 600; color: #374151; margin-bottom: 4px;">
                                {sample}
                            </div>
                            <div style="font-size: 20px; font-weight: 800; color: {gt_info['color']};">
                                {gt_info['icon']} {gt}
                            </div>
                        </div>
                    """, unsafe_allow_html=True)




def _display_grid_cards(filtered_dict, genotype_data):
    """Display genotypes in a responsive grid"""
    
    samples_list = list(filtered_dict.items())
    
    # Pre-build card HTML strings 
    card_htmls = []
    for sample, gt in samples_list:
        gt_info = genotype_data[gt]  
        card_htmls.append(f"""
            <div class="genotype-card" style="border-color: {gt_info['color']}; background: {gt_info['bg_color']};">
                <div class="card-sample-name">{sample}</div>
                <div class="card-genotype" style="color: {gt_info['color']};">{gt_info['icon']} {gt}</div>
                <div class="card-type">{gt_info['status'].replace('_', ' ').title()}</div>
            </div>
        """)
    
    # Build the complete HTML in one go
    grid_html = """
        <style>
            .genotype-grid {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(140px, 1fr));
                gap: 10px;
                margin: 14px 0;
            }
            .genotype-card {
                background: white;
                border: 2px solid #e5e7eb;
                border-radius: 10px;
                padding: 12px 8px;
                text-align: center;
                transition: all 0.25s ease;
                cursor: pointer;
                min-height: 100px;
                display: flex;
                flex-direction: column;
                justify-content: center;
            }
            .genotype-card:hover {
                transform: translateY(-2px) scale(1.03);
                box-shadow: 0 6px 20px rgba(0,0,0,0.15);
            }
            .card-sample-name {
                font-size: 13px;
                font-weight: 600;
                color: #374151;
                margin-bottom: 6px;
                word-break: break-word;
                line-height: 1.3;
            }
            .card-genotype {
                font-size: 22px;
                font-weight: 800;
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
                letter-spacing: 0.5px;
                margin-bottom: 4px;
                line-height: 1.2;
            }
            .card-type {
                font-size: 11px;
                color: #6B7280;
                text-transform: uppercase;
                letter-spacing: 0.5px;
                font-weight: 500;
            }
        </style>
        <div class="genotype-grid">
    """ + "".join(card_htmls) + "</div>"
    
    st.html(grid_html)

def display_motifs_as_bars_with_occurrences(sequence_name, motif_colors, motif_ids, spans, sequence, motif_names):
    # First calculate the motif occurrences
    motif_count = {}
    for motif_id in motif_ids:
        if motif_id != ".":
            motif_count[motif_id] = motif_count.get(motif_id, 0) + 1
    
    # Create two columns for side-by-side display
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Display the main bar visualization
        display_motifs_as_bars(sequence_name, motif_colors, motif_ids, spans, sequence, motif_names)
    
    with col2:
        if motif_count:
            st.markdown('<div class="occurrence-card">', unsafe_allow_html=True)
            
            # Create and display the occurrence chart
            bar_chart, df = plot_motif_bar(motif_count, motif_names, motif_colors, sequence_name)
            st.altair_chart(bar_chart, use_container_width=True)
            
            # Add summary statistics
            total_motifs = sum(motif_count.values())
            unique_motifs = len(motif_count)
            st.markdown(f"""
                <div class="occurrence-stats">
                    <div class="stat">
                        <div class="stat-value">{total_motifs}</div>
                        <div class="stat-label">Total Motifs</div>
                    </div>
                    <div class="stat">
                        <div class="stat-value">{unique_motifs}</div>
                        <div class="stat-label">Unique Types</div>
                    </div>
                </div>
            """, unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
    # Add the CSS for styling
    st.markdown("""
        <style>
            .occurrence-card {
                background: rgba(255, 255, 255, 0.95);
                backdrop-filter: blur(20px);
                border: 1px solid rgba(102, 126, 234, 0.15);
                border-radius: 16px;
                padding: 20px;
                margin: 20px 0;
                box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
            }
            .occurrence-stats {
                display: flex;
                gap: 15px;
                margin-top: 15px;
                justify-content: space-around;
            }
            .stat {
                text-align: center;
                background: #f8fafc;
                padding: 12px;
                border-radius: 12px;
                border: 1px solid #e2e8f0;
                flex: 1;
            }
            .stat-value {
                font-size: 18px;
                font-weight: 800;
                color: #667eea;
                margin-bottom: 4px;
            }
            .stat-label {
                font-size: 12px;
                color: #64748b;
                font-weight: 600;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }
        </style>
    """, unsafe_allow_html=True)

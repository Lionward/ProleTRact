import pandas as pd 
import re
import streamlit as st
import altair as alt

def parse_motif_range(motif_range):
    pattern = re.compile(r'\((\d+)-(\d+)\)')
    matches = pattern.findall(motif_range)
    ranges = [(int(start)-1, int(end)-1) for start, end in matches]
    return ranges



def motif_legend_html(motif_ids, motif_colors, motif_names):
    """
    Generate HTML for motif legend (excluding interruption legend)
    But now wraps it inside a tab labeled 'Motifs in region'
    """
    legend_html = ""
    unique_motifs = sorted(set(motif_ids))
    # Calculate motif sizes and prepare display names
    motif_display_data = []
    for motif_id in motif_colors.keys():
        color = motif_colors[int(motif_id)]
        motif_name = motif_names[int(motif_id)]
        motif_size = len(motif_name)
        motif_display_data.append({
            'id': motif_id,
            'color': color,
            'name': motif_name,
            'size': motif_size
        })
    
    # Sort by size to show larger motifs first
    motif_display_data.sort(key=lambda x: x['size'], reverse=True)
    
    for motif_data in motif_display_data:
        motif_name = motif_data['name']
        # Truncate very long motif names and show size
        if len(motif_name) > 20:
            display_name = f"{motif_name[:15]}... ({len(motif_name)}bp)"
        elif len(motif_name) > 10:
            display_name = f"{motif_name} ({len(motif_name)}bp)"
        else:
            display_name = f"{motif_name} ({len(motif_name)}bp)"
            
        legend_html += f"""
        <div class="legend-item" data-motif="{motif_data['id']}" title="{motif_name} - {len(motif_name)} bases">
            <span class="legend-color" style="background-color:{motif_data['color']};"></span>
            <span class="legend-text">{display_name}</span>
        </div>
        """
    
    legend_html += """
        <div class="legend-item" title="Interruption regions between motifs">
            <div class="legend-color interruption-color"></div>
            <span class="legend-text">Interruption</span>
        </div>
    """
    
    # Wrap legend in tab style like the sequence name, but tab labeled "Motifs in region"
    st.markdown("""
        <style>
        .motif-legend-container {
            margin: 16px 0;
            background: white;
            border-radius: 14px;
            box-shadow: 0 3px 14px rgba(0,0,0,0.08);
            overflow: hidden;
            border: 1px solid #e2e8f0;
        }
        .motif-legend-header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 10px 14px;
            font-weight: 700;
            font-size: 17px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .motif-count {
            background: rgba(255,255,255,0.18);
            padding: 3px 10px;
            border-radius: 14px;
            font-size: 15px;
            font-weight: 600;
        }
        .motif-legend-content {
            padding: 10px 10px 14px 10px;
            max-height: 135px;
            overflow-y: auto;
        }
        .motif-legend-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 7px;
        }
        .legend-item {
            display: flex;
            align-items: center;
            padding: 6px 9px;
            background: #f8fafc;
            border-radius: 8px;
            border: 2px solid transparent;
            transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
            cursor: pointer;
        }
        .legend-item:hover {
            transform: translateY(-1px);
            box-shadow: 0 4px 14px rgba(0,0,0,0.10);
            border-color: #667eea;
            background: white;
        }
        .legend-color {
            width: 15px;
            height: 15px;
            border-radius: 5px;
            margin-right: 7px;
            border: 2px solid white;
            box-shadow: 0 2px 6px rgba(0,0,0,0.09);
            flex-shrink: 0;
        }
        .interruption-color {
            background: linear-gradient(135deg, #fc8181, #e53e3e) !important;
        }
        .legend-text {
            font-size: 16px;
            font-weight: 600;
            color: #374151;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .motif-size-badge {
            background: #667eea;
            color: white;
            padding: 2px 6px;
            border-radius: 10px;
            font-size: 17px;
            font-weight: 700;
            margin-left: 7px;
        }
        .legend-stats {
            display: flex;
            gap: 8px;
            padding: 10px 10px 8px 10px;
            background: #f0f4ff;
            border-top: 1px solid #e2e8f0;
            font-size: 13px;
            color: #4b5563;
        }
        .stat-item {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .stat-value {
            background: white;
            padding: 2px 7px;
            border-radius: 7px;
            border: 1px solid #d1d5db;
            font-weight: 600;
            color: #1f2937;
        }
        /* Scrollbar styling */
        .motif-legend-content::-webkit-scrollbar {
            width: 5px;
        }
        .motif-legend-content::-webkit-scrollbar-thumb {
            background: #cbd5e0;
            border-radius: 7px;
        }
        .motif-legend-content::-webkit-scrollbar-track {
            background: #f1f5f9;
        }
        @media (max-width: 768px) {
            .motif-legend-grid {
                grid-template-columns: 1fr;
            }
        }
        </style>
    """, unsafe_allow_html=True)
    # Calculate statistics
    total_motifs = len(unique_motifs)
    motif_sizes = [len(motif_names[int(motif_id)]) for motif_id in unique_motifs]
    avg_size = sum(motif_sizes) / len(motif_sizes) if motif_sizes else 0
    max_size = max(motif_sizes) if motif_sizes else 0
    min_size = min(motif_sizes) if motif_sizes else 0
    
    # Create compact view for many motifs
    if total_motifs > 8:
        # Compact horizontal layout for many motifs
        compact_legend_html = ""
        for motif_data in motif_display_data[:12]:  # Show max 12 motifs
            motif_name = motif_data['name']
            display_name = f"{motif_name[:8]}..." if len(motif_name) > 8 else motif_name
            compact_legend_html += f"""
            <div class="compact-motif-item" data-motif="{motif_data['id']}" title="{motif_name} - {len(motif_name)} bases">
                <span class="compact-motif-color" style="background-color:{motif_data['color']};"></span>
                <span class="compact-motif-text">{display_name}</span>
            </div>
            """
        
        if total_motifs > 12:
            compact_legend_html += f"""
            <div class="compact-motif-item more-motifs">
                <span class="compact-motif-text">+{total_motifs - 12} more</span>
            </div>
            """

        
        st.html(f"""
            <style>
                .motif-legend-container {{
                    margin: 16px 0;
                    background: white;
                    border-radius: 12px;
                    box-shadow: 0 2px 12px rgba(0,0,0,0.06);
                    overflow: hidden;
                    border: 1px solid #e2e8f0;
                }}
                .motif-legend-header {{
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                    padding: 8px 12px;
                    font-weight: 700;
                    font-size: 14px;
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                }}
                .motif-count {{
                    background: rgba(255,255,255,0.2);
                    padding: 2px 8px;
                    border-radius: 12px;
                    font-size: 12px;
                    font-weight: 600;
                }}
                .compact-motifs-container {{
                    padding: 12px;
                    max-height: 80px;
                    overflow-y: auto;
                }}
                .compact-motifs-grid {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 6px;
                }}
                .compact-motif-item {{
                    display: flex;
                    align-items: center;
                    padding: 4px 8px;
                    background: #f8fafc;
                    border-radius: 6px;
                    border: 1px solid #e2e8f0;
                    transition: all 0.2s ease;
                    cursor: pointer;
                    font-size: 11px;
                }}
                .compact-motif-item:hover {{
                    transform: translateY(-1px);
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                    border-color: #667eea;
                }}
                .compact-motif-color {{
                    width: 10px;
                    height: 10px;
                    border-radius: 3px;
                    margin-right: 6px;
                    border: 1px solid white;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
                }}
                .compact-motif-text {{
                    font-weight: 600;
                    color: #374151;
                    white-space: nowrap;
                }}
                .more-motifs {{
                    background: #e5e7eb;
                    color: #6b7280;
                    font-style: italic;
                }}
                .legend-stats-compact {{
                    display: flex;
                    gap: 12px;
                    padding: 8px 12px;
                    background: #f0f4ff;
                    border-top: 1px solid #e2e8f0;
                    font-size: 11px;
                    color: #4b5563;
                }}
                .stat-item-compact {{
                    display: flex;
                    align-items: center;
                    gap: 4px;
                }}
                .stat-value-compact {{
                    background: white;
                    padding: 1px 6px;
                    border-radius: 4px;
                    border: 1px solid #d1d5db;
                    font-weight: 600;
                    color: #1f2937;
                }}
            </style>
            
            <div class="motif-legend-container">
                <div class="motif-legend-header">
                    <span>Motifs in Region</span>
                    <span class="motif-count">{total_motifs} motif{'' if total_motifs == 1 else 's'}</span>
                </div>
                
                <div class="compact-motifs-container">
                    <div class="compact-motifs-grid">
                        {compact_legend_html}
                    </div>
                </div>
                
                <div class="legend-stats-compact">
                    <div class="stat-item-compact">
                        <span>Range:</span>
                        <span class="stat-value-compact">{min_size}-{max_size}bp</span>
                    </div>
                    <div class="stat-item-compact">
                        <span>Avg:</span>
                        <span class="stat-value-compact">{avg_size:.1f}bp</span>
                    </div>
                </div>
            </div>
        """)
    else:
        # Original detailed view for few motifs
        st.html(f"""
            <div class="motif-legend-container">
                <div class="motif-legend-header">
                    <span>Motifs in Region</span>
                    <span class="motif-count">{total_motifs} motif{'' if total_motifs == 1 else 's'}</span>
                </div>
                
                <div class="legend-stats">
                    <div class="stat-item">
                        <span>Size Range:</span>
                        <span class="stat-value">{min_size}-{max_size}bp</span>
                    </div>
                    <div class="stat-item">
                        <span>Average:</span>
                        <span class="stat-value">{avg_size:.1f}bp</span>
                    </div>
                    <div class="stat-item">
                        <span>Total:</span>
                        <span class="stat-value">{total_motifs}</span>
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
    # Global override to upscale sequence visualization typography
    st.markdown("""
        <style>
            .sequence-dashboard, .sequence-dashboard * { font-size: 1.05rem !important; }
            .sequence-header { font-size: 1.08rem !important; }
            .sequence-length { font-size: 1.02rem !important; }
            .motif-legend-container, .motif-legend-container * { font-size: 1.05rem !important; }
            .motif-legend-header { font-size: 1.06rem !important; }
            .motif-count { font-size: 1.02rem !important; }
            .legend-stats .stat-item span, .legend-stats .stat-item .stat-value { font-size: 1.02rem !important; }
        </style>
    """, unsafe_allow_html=True)
    # Handle the situations where motifs are not available or sequence is just one base
    if (motif_ids == ["."]) or (isinstance(sequence, str) and len(sequence) <= 1):
        if sequence_name == "Ref":
            sequence_name += "seq"
        st.markdown(f"""
            <style>
                .sequence-dashboard {{
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
                    padding: 8px 12px;
                    font-weight: 700;
                    font-size: 12px;
                    display: flex;
                    gap: 10px;
                    align-items: center;
                }}
                .sequence-info-short {{
                    display: flex;
                    gap: 8px;
                    align-items: center;
                }}
                .sequence-length {{
                    font-size: 12px;
                    font-weight: 900;
                    color: #fffae3;
                    background: rgba(255,255,255,0.22);
                    padding: 1px 8px;
                    border-radius: 12px;
                    letter-spacing:0.5px;
                    box-shadow: 0 0px 2px rgba(0,0,0,0.04);
                }}
                .no-motif-section {{
                    padding: 12px;
                    background: linear-gradient(135deg, #fff7f7 0%, #feebeb 100%);
                    border: 1px solid #feb2b2;
                    border-radius: 8px;
                    margin: 12px;
                    text-align: center;
                }}
                .no-motif-icon {{
                    font-size: 17px;
                    margin-bottom: 6px;
                }}
                .no-motif-msg {{
                    color: #c53030;
                    font-weight: 600;
                    font-size: 16px;
                    margin-bottom: 6px;
                }}
                .no-motif-desc {{
                    color: #744210;
                    font-size: 10px;
                    opacity: 0.8;
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
                    height: 6px;
                }}
                .sequence-scroll-wrapper::-webkit-scrollbar-thumb {{
                    background: #cbd5e0;
                }}
                .sequence-content-wrapper {{
                    display: inline-block;
                    min-width: 100%;
                    overflow: visible;
                }}
                .sequence-content {{
                    padding: 10px;
                    white-space: nowrap;
                    overflow-x: visible;
                    overflow-y: hidden;
                    background: #f8fafc;
                    line-height: 1.6;
                    font-size: 14px;
                    font-weight: 500;
                }}
                .sequence-scale {{
                    position: relative;
                    height: 45px;
                    padding: 12px 10px 18px 10px;
                    font-size: 14px;
                    color: #718096;
                    font-weight: 600;
                    background: #f8fafc;
                    white-space: nowrap;
                    overflow: visible;
                    min-height: 45px;
                }}
                .sequence-scale .scale-marker {{
                    position: absolute;
                    transform: translateX(-50%);
                    white-space: nowrap;
                    z-index: 10;
                }}
                .sequence-scale .scale-marker:first-child {{ 
                    transform: none; 
                    left: 10px !important;
                }}
                .sequence-scale .scale-marker:last-child {{ 
                    transform: none; 
                    right: 10px !important;
                }}
                .sequence-content::-webkit-scrollbar {{
                    height: 6px;
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
                            <div class="sequence-content">
                                <span style="color:#4a5568; letter-spacing:0.5px;">{sequence}</span>
                            </div>
                            <div class="sequence-scale" id="sequence-scale-{sequence_name}">
                                <span class="scale-marker" style="left:10px;">0 bp</span>
                                <span class="scale-marker" style="left:25%;">{int(len(sequence) * 0.25)} bp</span>
                                <span class="scale-marker" style="left:50%;">{int(len(sequence) * 0.5)} bp</span>
                                <span class="scale-marker" style="left:75%;">{int(len(sequence) * 0.75)} bp</span>
                                <span class="scale-marker" style="right:10px;">{len(sequence)} bp</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        """, unsafe_allow_html=True)
        return

    # Process sequences with motifs
    ranges = parse_motif_range(spans)
    highlighted_sequence = ""
    previous_end = 0

    interruption_class = "interruption-segment motif-segment"
    interruption_style = "background: linear-gradient(135deg, #fc8181, #e53e3e); opacity: 0.85; border: 1px dashed rgba(255,255,255,0.5); color: #fff;"

    # Build highlighted sequence with interactive elements
    for idx, (start, end) in enumerate(ranges):
        motif = motif_ids[idx]
        color = motif_colors[int(motif)]
        motif_name = motif_names[int(motif)]

        # Add interruption if any
        if start > previous_end:
            interruption_sequence = sequence[previous_end:start]
            # Use the same outer box as motif, but with different color/border
            # Calculate 1-based positions for interruption
            int_start_pos_1based = previous_end + 1
            int_end_pos_1based = start
            highlighted_sequence += (
                f"<span class='{interruption_class}' data-type='interruption' "
                f"style='{interruption_style}' "
                f"data-content='⚡ <strong>Interruption Region</strong><br/>📍 Position: {int_start_pos_1based}-{int_end_pos_1based}<br/>📏 Length: {len(interruption_sequence)} bases<br/>🧪 Sequence: <code>{interruption_sequence}</code><br/>🔍 Type: Non-motif sequence'>"
                f"<span class='interruption-text'>{interruption_sequence}</span>"
                f"</span>"
            )

        # Add motif
        motif_sequence = sequence[start:end+1]
        motif_length = len(motif_sequence)
        
        # Calculate 1-based positions for display
        start_pos_1based = start + 1
        end_pos_1based = end + 1
        
        highlighted_sequence += (
            f"<span class='motif-segment motif-{motif} motif-length-{motif_length}' data-motif='{motif}' "
            f"style='background-color:{color};' "
            f"data-content='🧬 <strong>{motif_name}</strong><br/>📍 Position: {start_pos_1based}-{end_pos_1based}<br/>📏 Length: {motif_length} bases<br/>🧪 Sequence: <code>{motif_sequence}</code><br/>🆔 Motif ID: {motif}'>"
            f"<span class='motif-text'>{motif_sequence}</span>"
            f"</span>"
        )
        previous_end = end + 1

    # Add remaining sequence as interruption
    if previous_end < len(sequence):
        interruption_sequence = sequence[previous_end:]
        # Calculate 1-based positions for final interruption
        final_int_start_pos_1based = previous_end + 1
        final_int_end_pos_1based = len(sequence)
        highlighted_sequence += (
            f"<span class='{interruption_class}' data-type='interruption' "
            f"style='{interruption_style}' "
            f"data-content='⚡ <strong>Interruption Region</strong><br/>📍 Position: {final_int_start_pos_1based}-{final_int_end_pos_1based}<br/>📏 Length: {len(interruption_sequence)} bases<br/>🧪 Sequence: <code>{interruption_sequence}</code><br/>🔍 Type: Non-motif sequence'>"
            f"<span class='interruption-text'>{interruption_sequence}</span>"
            f"</span>"
        )

    if sequence_name == "Ref":
        sequence_name += "seq"

    # Calculate motif statistics
    legend_motif_sizes = [len(name) for name in motif_names]
    seen_sizes = set()
    ordered_unique_sizes = []
    for size in legend_motif_sizes:
        if size not in seen_sizes:
            ordered_unique_sizes.append(size)
            seen_sizes.add(size)
    motif_sizes_display = ", ".join([f"{size} bp" for size in ordered_unique_sizes])
    total_motifs = len(motif_ids)
    coverage = sum([end - start + 1 for start, end in ranges]) / len(sequence) * 100

    # Prepare supporting reads count html, show if supporting_reads is not None and is a number
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
            /* Keep all your existing CSS styles */
            .sequence-dashboard {{
                font-family: 'Inter', 'SF Pro Display', -apple-system, BlinkMacSystemFont, sans-serif;
                margin-bottom: 18px;
            }}
            .sequence-container {{
                border: 1px solid #e1e5e9;
                border-radius: 12px;
                overflow: hidden;
                box-shadow: 0 2px 10px rgba(0,0,0,0.05);
                background: white;
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
            .sequence-scroll-wrapper {{
                width: 100%;
                overflow-x: auto;
                overflow-y: visible;
                border-radius: 12px;
            }}
            .sequence-content-wrapper {{
                display: inline-block;
                min-width: 100%;
                overflow: visible;
            }}
            .sequence-content {{
                padding: 20px;
                white-space: nowrap;
                overflow-x: visible;
                overflow-y: hidden;
                background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%);
                line-height: 1.5;
                font-size: 13px;
                font-weight: 500;
                min-height: 42px;
            }}
            .sequence-scale {{
                position: relative;
                height: 45px;
                padding: 12px 20px 18px 20px;
                font-size: 14px;
                color: #718096;
                font-weight: 600;
                background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%);
                white-space: nowrap;
                overflow: visible;
                min-height: 45px;
            }}
            .sequence-scale .scale-marker {{
                position: absolute;
                transform: translateX(-50%);
                white-space: nowrap;
                z-index: 10;
            }}
            .sequence-scale .scale-marker:first-child {{ 
                transform: none; 
                left: 20px !important;
            }}
            .sequence-scale .scale-marker:last-child {{ 
                transform: none; 
                right: 20px !important;
            }}
            .motif-segment {{
                display: inline-block;
                padding: 3px 0px;
                margin: 0;
                border-radius: 10px;
                font-weight: 800;
                color: #ffffff;
                box-shadow: 0 4px 12px rgba(0,0,0,0.2);
                transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
                cursor: pointer;
                border: 3px solid rgba(255,255,255,0.5);
                position: relative;
                overflow: hidden;
                font-size: 14px;
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
                border-radius: 10px !important;
                margin: 0 !important;
                padding: 3px 0px !important;
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
                padding: 12px 16px;
                border-radius: 12px;
                font-size: 12px;
                pointer-events: none;
                opacity: 0;
                transform: translateY(15px) scale(0.95);
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                z-index: 1000;
                backdrop-filter: blur(20px);
                border: 2px solid rgba(255,255,255,0.2);
                max-width: 320px;
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
                padding: 2px 6px;
                border-radius: 4px;
                font-size: 11px;
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
                        <div class="sequence-content" id="sequence-content-{sequence_name}">
                            {highlighted_sequence}
                        </div>
                        <div class="sequence-scale" id="sequence-scale-{sequence_name}">
                            <span class="scale-marker" style="left:20px;">0 bp</span>
                            <span class="scale-marker" style="left:25%;">{int(len(sequence) * 0.25)} bp</span>
                            <span class="scale-marker" style="left:50%;">{int(len(sequence) * 0.5)} bp</span>
                            <span class="scale-marker" style="left:75%;">{int(len(sequence) * 0.75)} bp</span>
                            <span class="scale-marker" style="right:20px;">{len(sequence)} bp</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <script>
            // Simple global tooltip management
            let currentTooltip = null;
            
            function createTooltip() {{
                if (!currentTooltip) {{
                    currentTooltip = document.createElement('div');
                    currentTooltip.className = 'tooltip';
                    document.body.appendChild(currentTooltip);
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
            
            // Initialize tooltips for this sequence
            function initSequenceTooltips(sequenceId) {{
                const container = document.getElementById(sequenceId);
                if (!container) return;
                
                const segments = container.querySelectorAll('[data-content]');
                
                segments.forEach(segment => {{
                    // Mouse enter
                    segment.addEventListener('mouseenter', function(e) {{
                        const content = this.getAttribute('data-content');
                        showTooltip(content, e.pageX + 15, e.pageY + 15);
                    }});
                    
                    // Mouse leave
                    segment.addEventListener('mouseleave', hideTooltip);
                    
                    // Mouse move
                    segment.addEventListener('mousemove', function(e) {{
                        if (currentTooltip && currentTooltip.classList.contains('show')) {{
                            currentTooltip.style.left = (e.pageX + 15) + 'px';
                            currentTooltip.style.top = (e.pageY + 15) + 'px';
                        }}
                    }});
                }});
            }}
            
            // Sync scale width with sequence content width
            function syncScaleWidth() {{
                const sequenceContent = document.getElementById('sequence-content-{sequence_name}');
                const sequenceScale = document.getElementById('sequence-scale-{sequence_name}');
                if (sequenceContent && sequenceScale) {{
                    const contentWidth = sequenceContent.scrollWidth;
                    sequenceScale.style.width = contentWidth + 'px';
                    sequenceScale.style.minWidth = contentWidth + 'px';
                }}
            }}
            
            // Initialize when page loads
            document.addEventListener('DOMContentLoaded', function() {{
                initSequenceTooltips('sequence-content-{sequence_name}');
                syncScaleWidth();
            }});
            
            // Also try initializing after a short delay
            setTimeout(() => {{
                initSequenceTooltips('sequence-content-{sequence_name}');
                syncScaleWidth();
            }}, 100);
            
            // Sync on window resize
            window.addEventListener('resize', syncScaleWidth);
        </script>
    """)


def display_motifs_as_bars(sequence_name, motif_colors, motif_ids, spans, sequence, motif_names, supporting_reads=None):
    """
    This function draws motif bars for each motif in the sequence. 
    Motif information (name, count, coverage, sequence) is shown in a stylish tooltip on hover, not on the bar itself.
    This revision uses classic mouseover/mouseout events for tooltips and makes sure tooltips are attached per motif bar/interruption element.
    """
    from collections import Counter

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

    # Prepare supporting reads count html, show if supporting_reads is not None and is a number
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

    # Create interactive bar chart WITHOUT legend selection or visible legend
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

    # plot (fix label rotation on motif x-axis) and remove legend
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
    
    # Clean and normalize input
    gt_str = str(gt).strip()
    # For assembly mode, genotype may be a single allele ("0", "1") or a tuple/list
    if st.session_state.get("cohort_mode", "") == "assembly":
        if isinstance(gt, (list, tuple)):
            alleles = [str(a) for a in gt]
        else:
            # handle "0" or "1", or "0/1" (assembly data may sometimes use slash too)
            if "/" in gt_str:
                alleles = gt_str.split("/")
            elif "|" in gt_str:
                alleles = gt_str.split("|")
            else:
                alleles = [gt_str]
    else:
        # Default - try splitting with slash or pipe; else treat as single allele
        if "/" in gt_str:
            alleles = gt_str.split("/")
        elif "|" in gt_str:
            alleles = gt_str.split("|")
        else:
            alleles = [gt_str]

    # Remove empty alleles (shouldn't happen with good VCFs)
    alleles = [a for a in alleles if a not in [".", ""]]

    # Accept single-haplotype genotypes as well (chrX, chrY, assemblies)
    if len(alleles) == 1:
        allele = alleles[0]
        if allele == "0":
            return {
                'description': 'Hemizygous Reference (0)',
                'interpretation': 'Only one haplotype; matches reference',
                'color': '#10B981',
                'bg_color': '#D1FAE5',
                'icon': '🟢',
                'status': 'hemizygous_ref'
            }
        else:
            return {
                'description': f'Hemizygous Alternative ({allele})',
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
                'description': 'Homozygous Reference (0/0)',
                'interpretation': 'Both haplotypes identical to reference',
                'color': '#10B981',
                'bg_color': '#D1FAE5',
                'icon': '🟢',
                'status': 'homozygous_ref'
            }
        else:
            return {
                'description': f'Homozygous Alternative ({allele1}/{allele2})',
                'interpretation': f'Both haplotypes different from reference (allele {allele1})',
                'color': '#F59E0B',
                'bg_color': '#FEF3C7',
                'icon': '🟡',
                'status': 'homozygous_alt'
            }
    else:
        if allele1 == "0" or allele2 == "0":
            return {
                'description': f'Heterozygous Reference/Alternative ({allele1}/{allele2})',
                'interpretation': 'One haplotype like reference, one different',
                'color': '#3B82F6',
                'bg_color': '#DBEAFE',
                'icon': '🔵',
                'status': 'heterozygous_ref_alt'
            }
        else:
            return {
                'description': f'Heterozygous Alternative ({allele1}/{allele2})',
                'interpretation': f'Both haplotypes different from reference (alleles {allele1} and {allele2})',
                'color': '#8B5CF6',
                'bg_color': '#EDE9FE',
                'icon': '🟣',
                'status': 'heterozygous_alt_alt'
            }


def display_genotype_card(gt, sample_name="Sample", show_details=True):
    """
    Display genotype information in a creative card format
    
    Args:
        gt (str): Genotype string
        sample_name (str): Name of the sample
        show_details (bool): Whether to show detailed interpretation
    """
    interpretation = interpret_genotype(gt)
    
    st.html(f"""
        <style>
            .genotype-card {{
                background: linear-gradient(135deg, {interpretation['bg_color']} 0%, rgba(255,255,255,0.9) 100%);
                border: 1px solid {interpretation['color']};
                border-radius: 8px;
                padding: 8px 12px;
                margin: 8px 0;
                box-shadow: 0 2px 8px rgba(0,0,0,0.06);
                transition: all 0.2s ease;
                position: relative;
                overflow: hidden;
            }}
            .genotype-card:hover {{
                transform: translateY(-1px);
                box-shadow: 0 3px 12px rgba(0,0,0,0.1);
            }}
            .genotype-content {{
                display: flex;
                align-items: center;
                gap: 8px;
                flex-wrap: wrap;
            }}
            .genotype-icon {{
                font-size: 20px;
            }}
            .genotype-title {{
                font-size: 20px;
                font-weight: 700;
                color: {interpretation['color']};
                margin: 0;
            }}
            .genotype-value {{
                background: {interpretation['color']};
                color: white;
                padding: 2px 6px;
                border-radius: 8px;
                font-size: 20px;
                font-weight: 800;
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
                letter-spacing: 0.3px;
            }}
            .genotype-description {{
                font-size: 20px;
                font-weight: 600;
                color: #374151;
                margin: 0;
            }}
            .genotype-stats {{
                display: flex;
                gap: 20px;
                margin-left: auto;
            }}
            .stat-item {{
                font-size: 20px;
                color: #6B7280;
                font-weight: 500;
            }}
            .stat-value {{
                color: {interpretation['color']};
                font-weight: 700;
            }}
        </style>
        
        <div class="genotype-card">
            <div class="genotype-content">
                <span class="genotype-icon">{interpretation['icon']}</span>
                <span class="genotype-title">{sample_name}:</span>
                <span class="genotype-value">{gt}</span>
                <span class="genotype-description">{interpretation['description']}</span>
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


def create_genotype_comparison_matrix(genotypes_dict):
    """
    Create a visual comparison matrix of genotypes across samples, hidden behind an expander
    Args:
        genotypes_dict (dict): Dictionary with sample names as keys and genotypes as values
    """


    st.markdown("### 🧬 Genotype Comparison Matrix")

    # Create comparison data
    samples = list(genotypes_dict.keys())
    unique_genotypes = list(set(genotypes_dict.values()))

    # Color mapping for different genotypes
    genotype_colors = {}
    for i, gt in enumerate(unique_genotypes):
        interpretation = interpret_genotype(gt)
        genotype_colors[gt] = interpretation['color']

    # Adjusted/Smaller matrix CSS
    matrix_html = """
        <style>
            .genotype-matrix {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
                gap: 8px;
                margin: 14px 0;
            }
            .matrix-cell {
                background: white;
                border: 1.5px solid #e5e7eb;
                border-radius: 8px;
                padding: 7px 5px 9px 5px;
                text-align: center;
                transition: all 0.22s ease;
                position: relative;
                overflow: hidden;
                min-width: 0;
                min-height: 0;
            }
            .matrix-cell:hover {
                transform: translateY(-1px) scale(1.02);
                box-shadow: 0 4px 12px rgba(0,0,0,0.10);
            }
            .sample-name {
                font-size: 14px;
                font-weight: 600;
                color: #374151;
                margin-bottom: 3px;
                white-space: pre-line;
                overflow-wrap: anywhere;
                word-break: break-all;
            }
            .genotype-display {
                font-size: 18px;
                font-weight: 800;
                font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
                letter-spacing: 0.6px;
                margin-bottom: 2px;
                line-height: 1.1;
            }
            .genotype-type {
                font-size: 16px;
                color: #6B7280;
                text-transform: uppercase;
                letter-spacing: 0.28px;
            }
        </style>
        <div class="genotype-matrix">
    """

    for sample, gt in genotypes_dict.items():
        interpretation = interpret_genotype(gt)
        matrix_html += f"""
            <div class="matrix-cell" style="border-color: {interpretation['color']}; background: {interpretation['bg_color']};">
                <div class="sample-name">{sample}</div>
                <div class="genotype-display" style="color: {interpretation['color']};">{interpretation['icon']} {gt}</div>
                <div class="genotype-type">{interpretation['status'].replace('_', ' ').title()}</div>
            </div>
        """

    matrix_html += "</div>"

    st.html(matrix_html)

    # Add summary statistics
    st.markdown("#### 📊 Genotype Summary")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Samples", len(samples))

    with col2:
        st.metric("Unique Genotypes", len(unique_genotypes))

    with col3:
        homozygous_count = sum(1 for gt in genotypes_dict.values() if interpret_genotype(gt)['status'].startswith('homozygous'))
        st.metric("Homozygous", homozygous_count)

    with col4:
        heterozygous_count = sum(1 for gt in genotypes_dict.values() if interpret_genotype(gt)['status'].startswith('heterozygous'))
        st.metric("Heterozygous", heterozygous_count)

def display_motifs_as_bars_with_occurrences(sequence_name, motif_colors, motif_ids, spans, sequence, motif_names):
    # First calculate motif occurrences
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
    # Add this CSS for styling
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

"""
Utility functions for visualization and data processing.
"""
import streamlit as st
import math
import colorsys
from . import parsers


def get_color_palette(num_colors):
    
    # Primary palette matching your app's gradient theme
    primary_palette = [
        '#667eea',  # Primary blue
        '#764ba2',  # Primary purple
        '#4fd1c7',  # Teal
        '#68d391',  # Green
        '#f6ad55',  # Orange
        '#fc8181',  # Coral
        '#f093fb',  # Pink
        '#f6e05e',  # Yellow
    ]
    
    # Secondary palette - darker shades
    secondary_palette = [
        '#5a67d8', '#6b46c1', '#38a169', '#dd6b20',
        '#e53e3e', '#d53f8c', '#d69e2e', '#319795'
    ]
    
    # Tertiary palette - lighter shades
    tertiary_palette = [
        '#a3bffa', '#b794f4', '#9ae6b4', '#fbd38d',
        '#fed7d7', '#fbb6ce', '#faf089', '#81e6d9'
    ]
    
    # Combine all palettes
    combined_palette = primary_palette + secondary_palette + tertiary_palette
    
    # If we need more colors, generate them algorithmically
    if num_colors > len(combined_palette):
        generated_colors = []
        # Use golden ratio for optimal color distribution
        golden_ratio_conjugate = (math.sqrt(5) - 1) / 2
        
        for i in range(num_colors - len(combined_palette)):
            # Generate colors with good separation in HSV space
            hue = (i * golden_ratio_conjugate) % 1.0
            # Vary saturation and value for better distinction
            saturation = 0.6 + 0.2 * (i % 3) / 3
            value = 0.7 + 0.2 * ((i + 1) % 2)
            
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            hex_color = '#{:02x}{:02x}{:02x}'.format(
                int(rgb[0] * 255), 
                int(rgb[1] * 255), 
                int(rgb[2] * 255)
            )
            generated_colors.append(hex_color)
        
        return combined_palette + generated_colors[:num_colors]
    
    return combined_palette[:num_colors]


def get_results_cohort(region, files, file_paths, cohort_mode, parse_record_func, parse_record_assembly_func):
    """
    Get cohort results by parsing records from files.
    
    Args:
        region: Genomic region string
        files: List of VCF files
        file_paths: List of file paths
        cohort_mode: Mode ('assembly' or 'reads')
        parse_record_func: Function to parse reads-based records
        parse_record_assembly_func: Function to parse assembly records
    
    Returns:
        dict: Dictionary of sample names to records
    """
    samples_results = {}
    chrom, start_end = region.split(":")
    start, end = start_end.split("-")
    start = int(start) - 1
    end = int(end) - 1
    region = f"{chrom}:{start}-{end}"

    for i in range(len(files)):
        sample_name = file_paths[i].split(".")[0]
        if cohort_mode == "assembly":
            record = parse_record_assembly_func(files[i], region)
        else:
            record = parse_record_func(files[i], region)
        samples_results[sample_name] = record
    return samples_results


def get_results_hgsvc_pop(region, files, file_paths, parse_record_assembly_func):
    """
    Get HGSVC population results.
    
    Args:
        region: Genomic region string
        files: List of VCF files
        file_paths: List of file paths
        parse_record_assembly_func: Function to parse assembly records
    
    Returns:
        dict: Dictionary of sample names to records, or None if files not available
    """
    if files is None:
        return None
    samples_results = {}
    for i in range(len(files)):
        sample_name = file_paths[i].split(".")[0]
        record = parse_record_assembly_func(files[i], region)
        samples_results[sample_name] = record
    return samples_results


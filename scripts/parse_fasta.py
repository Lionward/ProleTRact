import pysam
import pandas as pd 
# import numpy as np

def parse_fasta_reads_cut(fasta_file, seq_col_name="seq"):
    # Initialize lists to store data
    chr_list = []
    start_list = []
    end_list = []
    seq_list = []
    
    # Open the FASTA file using pysam
    fasta = pysam.FastaFile(fasta_file)
    
    # Iterate over each reference (read name) in the FASTA file
    for read_name in fasta.references:
        # Extract the sequence for this read
        seq = fasta.fetch(read_name)
        
        # Find the portion before the underscore (_) to extract chr, start, and end
        main_part = read_name.split('_')[0]
        
        # Parse the main part as chr:start-end
        chr_part, region = main_part.split(':')
        start, end = region.split('-')
        
        # Append chr, start, end, and seq to the respective lists
        chr_list.append(chr_part)
        start_list.append(int(start))
        end_list.append(int(end))
        seq_list.append(seq)
    
    # Close the FASTA file
    fasta.close()
    
    # Create a DataFrame
    df = pd.DataFrame({
        'chrom': chr_list,
        'start': start_list,
        'end': end_list,
        'seq': seq_list
    })

    df = df.rename(columns={'seq': seq_col_name})
    
    return df


# Function to count non-overlapping occurrences of each motif in a sequence
def count_motifs(sequence, motifs):
    counts = []
    for motif in motifs:
        count = sequence.count(motif)  # Count non-overlapping occurrences
        counts.append(count)
    return tuple(counts)


def compare_with_asm(tr_df, asm_df):
    # Step 1: Select the desired columns from tr_df
    selected_columns = tr_df[['chrom', 'start', 'end', 'motifs', 'motif_count_h1', 'motif_count_h2']]

    # Step 2: Merge the selected DataFrame with asm_df
    merged_df = pd.merge(asm_df, selected_columns, on=['chrom', 'start', 'end'], how='inner')

    # Apply the count_motifs function for asm_h1_seq and asm_h2_seq
    merged_df['asm_h1_motif_count'] = merged_df.apply(lambda row: count_motifs(row['asm_h1_seq'], row['motifs']), axis=1)
    merged_df['asm_h2_motif_count'] = merged_df.apply(lambda row: count_motifs(row['asm_h2_seq'], row['motifs']), axis=1)

    return merged_df
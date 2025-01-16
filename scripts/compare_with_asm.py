import pandas as pd
import numpy as np
import Levenshtein
from pandarallel import pandarallel

# Define functions for similarity calculations
def normalized_levenshtein(seq1, seq2):
    """Calculate Normalized Levenshtein Similarity."""
    edit_dist = Levenshtein.distance(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    return 1 - (edit_dist / max_len)

def relative_length_difference(seq1, seq2):
    """Calculate Relative Length Difference."""
    len1, len2 = len(seq1), len(seq2)
    return 1 - abs(len1 - len2) / max(len1, len2)

# Example DataFrame
data = {
    'motif_concat_h1': ["AACCTGCTGCTTCCTGGAGGAAGACAGTCCCTCTGTCCCTCTGTCTA"],
    'motif_concat_h2': ["AACCTGCTGCTTCCTGGAGGAAGACAGTCCCTCTGTCCCTCTGTCT"],
    'asm_h1_seq': ["AACCTGCTGCTTCCTGGAGGAAGACAGTCCCTCTGTCCCTCTGTCTA"],
    'asm_h2_seq': ["AACCTGCTGCTTCCTGGAGGAAGACAGTCCCTCTGTCCCTCTGTCT"]
}

merged_df = pd.DataFrame(data)


# print(merged_df)
# print(merged_df[['seq_sim_11',  'len_sim_11',  'seq_sim_12' , 'len_sim_12', 'seq_sim_21',  'len_sim_21',  'seq_sim_22' , 'len_sim_22', 'seq_sim', 'len_sim', 'assigned_haplotype']])

def compare_sample_to_asm(df, asm_df):
    # merge the two dataframes
    merged_df = pd.merge(df, asm_df, on=['chrom', 'start', 'end'], how='inner')

    # Preallocate columns for results
    new_cols = ['seq_sim_11', 'len_sim_11', 'seq_sim_12', 'len_sim_12',
                'seq_sim_21', 'len_sim_21', 'seq_sim_22', 'len_sim_22',
                'seq_sim', 'len_sim', 'assigned_haplotype']

    for col in new_cols:
        if col == 'assigned_haplotype':
            merged_df[col] = None  # Set dtype to object for string compatibility
        else:
            merged_df[col] = np.nan

    # Compute similarities row-wise
    for index, row in merged_df.iterrows():
        seq1_h1 = row['motif_concat_h1']
        seq2_h1 = row['asm_h1_seq']
        seq2_h2 = row['asm_h2_seq']
        seq1_h2 = row['motif_concat_h2']

        # Calculate similarities
        seq_sim_11 = normalized_levenshtein(seq1_h1, seq2_h1)
        len_sim_11 = relative_length_difference(seq1_h1, seq2_h1)

        seq_sim_12 = normalized_levenshtein(seq1_h1, seq2_h2)
        len_sim_12 = relative_length_difference(seq1_h1, seq2_h2)

        seq_sim_21 = normalized_levenshtein(seq1_h2, seq2_h1)
        len_sim_21 = relative_length_difference(seq1_h2, seq2_h1)

        seq_sim_22 = normalized_levenshtein(seq1_h2, seq2_h2)
        len_sim_22 = relative_length_difference(seq1_h2, seq2_h2)

        merged_df.loc[index, 'seq_sim_11'] = seq_sim_11
        merged_df.loc[index, 'len_sim_11'] = len_sim_11

        merged_df.loc[index, 'seq_sim_12'] = seq_sim_12
        merged_df.loc[index, 'len_sim_12'] = len_sim_12

        merged_df.loc[index, 'seq_sim_21'] = seq_sim_21
        merged_df.loc[index, 'len_sim_21'] = len_sim_21

        merged_df.loc[index, 'seq_sim_22'] = seq_sim_22
        merged_df.loc[index, 'len_sim_22'] = len_sim_22

        # Determine the best match
        if seq_sim_11 + seq_sim_22 > seq_sim_12 + seq_sim_21:
            merged_df.loc[index, 'seq_sim'] = (seq_sim_11 + seq_sim_22) / 2
            merged_df.loc[index, 'len_sim'] = (len_sim_11 + len_sim_22) / 2
            merged_df.loc[index, 'assigned_haplotype'] = "11_22"
        else:
            merged_df.loc[index, 'seq_sim'] = (seq_sim_12 + seq_sim_21) / 2
            merged_df.loc[index, 'len_sim'] = (len_sim_12 + len_sim_21) / 2
            merged_df.loc[index, 'assigned_haplotype'] = "12_21"
        
    return merged_df



# Initialize pandarallel with progress bars
pandarallel.initialize(progress_bar=True, nb_workers = 4)

def compare_sample_to_asm_parallel(df, asm_df):
    # Merge the two dataframes
    merged_df = pd.merge(df, asm_df, on=['chrom', 'start', 'end'], how='inner')

    # Preallocate columns for results
    new_cols = ['seq_sim_11', 'len_sim_11', 'seq_sim_12', 'len_sim_12',
                'seq_sim_21', 'len_sim_21', 'seq_sim_22', 'len_sim_22',
                'seq_sim', 'len_sim', 'assigned_haplotype']

    for col in new_cols:
        if col == 'assigned_haplotype':
            merged_df[col] = None  # Set dtype to object for string compatibility
        else:
            merged_df[col] = np.nan
    
    def compute_similarities(row):
        seq1_h1 = row['motif_concat_h1']
        seq2_h1 = row['asm_h1_seq']
        seq2_h2 = row['asm_h2_seq']
        seq1_h2 = row['motif_concat_h2']
        
        # Calculate similarities
        results = {}
        results['seq_sim_11'] = normalized_levenshtein(seq1_h1, seq2_h1)
        results['len_sim_11'] = relative_length_difference(seq1_h1, seq2_h1)
        
        results['seq_sim_12'] = normalized_levenshtein(seq1_h1, seq2_h2)
        results['len_sim_12'] = relative_length_difference(seq1_h1, seq2_h2)
        
        results['seq_sim_21'] = normalized_levenshtein(seq1_h2, seq2_h1)
        results['len_sim_21'] = relative_length_difference(seq1_h2, seq2_h1)
        
        results['seq_sim_22'] = normalized_levenshtein(seq1_h2, seq2_h2)
        results['len_sim_22'] = relative_length_difference(seq1_h2, seq2_h2)
        
        if results['seq_sim_11'] + results['seq_sim_22'] > results['seq_sim_12'] + results['seq_sim_21']:
            results['seq_sim'] = (results['seq_sim_11'] + results['seq_sim_22']) / 2
            results['len_sim'] = (results['len_sim_11'] + results['len_sim_22']) / 2
            results['assigned_haplotype'] = "11_22"
        else:
            results['seq_sim'] = (results['seq_sim_12'] + results['seq_sim_21']) / 2
            results['len_sim'] = (results['len_sim_12'] + results['len_sim_21']) / 2
            results['assigned_haplotype'] = "12_21"
        
        return pd.Series(results)

    # Use pandarallel to apply computations in parallel
    result_df = merged_df.parallel_apply(compute_similarities, axis=1)
    merged_df[new_cols] = result_df
    
    return merged_df
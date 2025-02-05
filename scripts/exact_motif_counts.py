def filter_ambiguous_motifs_with_indices(motifs):
    # Set to keep track of ambiguous indices
    # two motifs are ambiguous if one is a substr of the other
    ambiguous_indices = set()
    
    for i, motif in enumerate(motifs):
        for j, other_motif in enumerate(motifs):
            if i != j and motif in other_motif:
                ambiguous_indices.add(i)
                ambiguous_indices.add(j)

    # Get the indices of non-ambiguous motifs
    non_ambiguous_indices = [i for i in range(len(motifs)) if i not in ambiguous_indices]
    
    return non_ambiguous_indices

# Example usage
# motifs = ['CTTTTTT', 'CTTTTT', 'TTTTT', 'AAAAAA', 'AA']
# filtered_indices = filter_ambiguous_motifs_with_indices(motifs)
# print(filtered_indices)  # Output will be indices of non-ambiguous motifs


def mean_abs_err(l1, l2, indices):
    # print(f'l1 = {l1}, l2 = {l2}')
    # Ensure both inputs are lists or tuples
    if isinstance(l1, (list, tuple)) and isinstance(l2, (list, tuple)):
        # n = len(l1)
        err = [abs(l1[i] - l2[i]) for i in indices]
        return sum(err) / len(indices)
    else:
        return np.nan  # Return NaN if inputs are not valid


def compare_counts_with_asm(row):
    # returns the dist with asm counts, and if the input motifs are ambiguous (all)
    non_ambig_motifs = filter_ambiguous_motifs_with_indices(row['motifs'])

    if not non_ambig_motifs:
        # non ambigous motif set is empty
        return np.nan, np.nan, True

    # Calculate mean absolute errors
    d11 = mean_abs_err(row['motif_count_h1'], row['asm_h1_motif_count'], non_ambig_motifs)
    d22 = mean_abs_err(row['motif_count_h2'], row['asm_h2_motif_count'], non_ambig_motifs)

    d12 = mean_abs_err(row['motif_count_h1'], row['asm_h2_motif_count'], non_ambig_motifs)
    d21 = mean_abs_err(row['motif_count_h2'], row['asm_h1_motif_count'], non_ambig_motifs)

    if np.isnan(d11) or np.isnan(d22) or np.isnan(d12) or np.isnan(d21):
        return np.nan, np.nan, False

    total_asm_motif_count = sum([row['asm_h1_motif_count'][i] for i in non_ambig_motifs]) + sum([row['asm_h2_motif_count'][i] for i in non_ambig_motifs])

    abs_MAEs = min(d11 + d22, d12 + d21)
    relative_MAEs = abs_MAEs if total_asm_motif_count == 0 else abs_MAEs / total_asm_motif_count
    
    return abs_MAEs, relative_MAEs, False
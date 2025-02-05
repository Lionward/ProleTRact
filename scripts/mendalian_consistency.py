########## Calulate the mendalian consistency

def is_mendelian_consistent_chrY(row, tolerance=0):
    '''
    checks if the child's genotype on chromosome Y is consistent with the father's genotype, allowing for an off-by-one error
    
    Args:
        row: pandas series
        tolerance: int, the allowed deviation (default is 0)
    Returns:

        bool
    '''

    def is_close(a, b):
    
        return abs(a - b) <= tolerance

    # Only males have a Y chromosome
    if row['sex_child'] == 'male':
        child_h1 = row[f'CN_H1_child']
        father_h1 = row[f'CN_H1_father']
        
        return is_close(child_h1, father_h1)
    
    # Females do not have a Y chromosome
    return False

def is_mendelian_consistent_chrX(row, tolerance=0):
    '''
    checks if the child's genotype on chromosome X is consistent with the parents' genotypes, allowing for an off-by-one error
    
    Args:
        row: pandas series
        tolerance: int, the allowed deviation (default is 0)
    Returns:
        bool
    '''

    child_h1 = row[f'CN_H1_child']
    child_h2 = row[f'CN_H2_child']
    father_h1 = row[f'CN_H1_father']
    mother_h1 = row[f'CN_H1_mother']
    mother_h2 = row[f'CN_H2_mother']
    
    def is_close(a, b):
        return abs(a - b) <= tolerance
    
    # If the child is male (XY), he inherits his X chromosome from his mother
    if row['sex_child'] == 'male':
        return is_close(child_h1, mother_h1) or is_close(child_h1, mother_h2)
    
    # If the child is female (XX), she inherits one X chromosome from each parent
    combinations = [
        (is_close(child_h1, father_h1) and (is_close(child_h2, mother_h1) or is_close(child_h2, mother_h2))),
        (is_close(child_h2, father_h1) and (is_close(child_h1, mother_h1) or is_close(child_h1, mother_h2)))
    ]

    return any(combinations)

def is_mendelian_consistent(row, tolerance=0):
    '''
    checks if the child's genotype is consistent with the parents' genotypes, allowing for an off-by-one error
    
    Args:
        row: pandas series
    Returns:
        bool
    '''
    if row['chrom'] == 'chrX':
        return is_mendelian_consistent_chrX(row, tolerance)
    if row['chrom'] == 'chrY':
        return is_mendelian_consistent_chrY(row, tolerance)
    
    child_h1 = row[f'CN_H1_child']
    child_h2 = row[f'CN_H2_child']
    father_h1 = row[f'CN_H1_father']
    father_h2 = row[f'CN_H2_father']
    mother_h1 = row[f'CN_H1_mother']
    mother_h2 = row[f'CN_H2_mother']
    
    def is_close(a, b):
        return abs(a - b) <= tolerance
    
    combinations = [
        (is_close(child_h1, father_h1) and is_close(child_h2, mother_h1)),
        (is_close(child_h1, father_h1) and is_close(child_h2, mother_h2)),
        (is_close(child_h1, father_h2) and is_close(child_h2, mother_h1)),
        (is_close(child_h1, father_h2) and is_close(child_h2, mother_h2)),
        (is_close(child_h2, father_h1) and is_close(child_h1, mother_h1)),
        (is_close(child_h2, father_h1) and is_close(child_h1, mother_h2)),
        (is_close(child_h2, father_h2) and is_close(child_h1, mother_h1)),
        (is_close(child_h2, father_h2) and is_close(child_h1, mother_h2))
    ]

    return any(combinations)

def count_motifs(motif_ids, motifs):
    # Efficiently count motifs using numpy and return a count dictionary
    if pd.isna(motif_ids) or pd.isna(motifs):
        return {}
    unique, counts = np.unique(motif_ids, return_counts=True)
    return dict(zip([motifs[int(i)] for i in unique], counts))


def is_mendelian_consistent_based_on_motif_chrY(row):
    '''
    checks if the child's genotype on chromosome Y is consistent with the father's genotype
    
    Args:
        row: pandas series
    Returns:
        bool
    '''
    if row['sex_child'] == 'male':
        motifs_child = row['motifs_child'] if not pd.isna(row['motifs_child']) else []
        motifs_father = row['motifs_father'] if not pd.isna(row['motifs_father']) else []
        child_h1 = row[f'motif_ids_H1_child']
        child_h1 = {motifs_child[int(i)]: child_h1.count(str(i)) for i in child_h1}
        father_h1 = row[f'motif_ids_H1_father']
        father_h1 = {motifs_father[int(i)]: father_h1.count(str(i)) for i in father_h1}
        


        return child_h1 == father_h1 
    
    return False

def is_mendelian_consistent_based_on_motif_chrX(row):
    '''
    checks if the child's genotype on chromosome X is consistent with the parents' genotypes
    
    Args:
        row: pandas series
    Returns:
        bool
    '''
    motifs_child=  row['motifs_child'] if not pd.isna(row['motifs_child']) else []
    motifs_father = row['motifs_father'] if not pd.isna(row['motifs_father']) else []
    motifs_mother = row['motifs_mother'] if not pd.isna(row['motifs_mother']) else []
    
    child_h1 = row[f'motif_ids_H1_child']
    child_h1 = {motifs_child[int(i)]: child_h1.count(str(i)) for i in child_h1}
    child_h2 = row[f'motif_ids_H2_child']
    child_h2 = {motifs_child[int(i)]: child_h2.count(str(i)) for i in child_h2}
    father_h1 = row[f'motif_ids_H1_father']
    father_h1 = {motifs_father[int(i)]: father_h1.count(str(i)) for i in father_h1}
    mother_h1 = row[f'motif_ids_H1_mother']
    mother_h1 = {motifs_mother[int(i)]: mother_h1.count(str(i)) for i in mother_h1}
    mother_h2 = row[f'motif_ids_H2_mother']
    mother_h2 = {motifs_mother[int(i)]: mother_h2.count(str(i)) for i in mother_h2}

    if row['sex_child'] == 'male':
        return child_h1 in [mother_h1, mother_h2] or child_h2 in [mother_h1, mother_h2]

    combinations = [
        (child_h1 == father_h1 and child_h2 in [mother_h1, mother_h2]),
        (child_h2 == father_h1 and child_h1 in [mother_h1, mother_h2])
    ]
    return any(combinations)



def is_mendelian_consistent_based_on_motif(row):
    '''
    checks if the child's genotype is consistent with the parents' genotypes
    
    Args:
        row: pandas series
    Returns:
        bool
    '''
    if row['chrom'] == 'chrX':
        return is_mendelian_consistent_based_on_motif_chrX(row)
    if row['chrom'] == 'chrY':
        return is_mendelian_consistent_based_on_motif_chrY(row)
    motifs_child=  row['motifs_child'] if not pd.isna(row['motifs_child']) else []
    motifs_father = row['motifs_father'] if not pd.isna(row['motifs_father']) else []
    motifs_mother = row['motifs_mother'] if not pd.isna(row['motifs_mother']) else []
    
    child_h1 = row[f'motif_ids_H1_child']
    child_h1 = {motifs_child[int(i)]: child_h1.count(str(i)) for i in child_h1}
    child_h2 = row[f'motif_ids_H2_child']
    child_h2 = {motifs_child[int(i)]: child_h2.count(str(i)) for i in child_h2}
    father_h1 = row[f'motif_ids_H1_father']
    father_h1 = {motifs_father[int(i)]: father_h1.count(str(i)) for i in father_h1}
    father_h2 = row[f'motif_ids_H2_father']
    father_h2 = {motifs_father[int(i)]: father_h2.count(str(i)) for i in father_h2}
    mother_h1 = row[f'motif_ids_H1_mother']
    mother_h1 = {motifs_mother[int(i)]: mother_h1.count(str(i)) for i in mother_h1}
    mother_h2 = row[f'motif_ids_H2_mother']
    mother_h2 = {motifs_mother[int(i)]: mother_h2.count(str(i)) for i in mother_h2}
    
    
    combinations = [
        (child_h1 == father_h1 and child_h2 == mother_h1),
        (child_h1 == father_h1 and child_h2 == mother_h2),
        (child_h1 == father_h2 and child_h2 == mother_h1),
        (child_h1 == father_h2 and child_h2 == mother_h2),
        (child_h2 == father_h1 and child_h1 == mother_h1),
        (child_h2 == father_h1 and child_h1 == mother_h2),
        (child_h2 == father_h2 and child_h1 == mother_h1),
        (child_h2 == father_h2 and child_h1 == mother_h2)
    ]

    return any(combinations)


def get_mendelian_consistency(child_df, father_df, mother_df, tool, child_sex):
    '''
    This function merges the child, father and mother dataframes and calculates the mendelian consistency for each region
    
    Args:
        child_df: pandas dataframe
        father_df: pandas dataframe
        mother_df: pandas dataframe
        tool: str, the tool used to call the variants
    Returns:
        merged_df: pandas dataframe

    '''
    merged_df_1 = child_df.merge(father_df, on=['chrom', 'start', 'end'], how='left', suffixes=('_child', '_father'))
    merged_df_2 = child_df.merge(mother_df, on=['chrom', 'start', 'end'], how='left', suffixes=('_child', '_mother'))
    merged_df = merged_df_1.merge(merged_df_2, on=['chrom', 'start', 'end', 'motifs_child', 'motif_ids_H1_child', 'motif_ids_H2_child', 'CN_H1_child', 'CN_H2_child', 'CN_ref_child', 'GT_child'], how='left')
    # replace the nan valuse in the cn columns with 0
    merged_df['CN_H1_father'] = merged_df['CN_H1_father'].fillna(0)
    merged_df['CN_H2_father'] = merged_df['CN_H2_father'].fillna(0)
    merged_df['CN_H1_mother'] = merged_df['CN_H1_mother'].fillna(0)
    merged_df['CN_H2_mother'] = merged_df['CN_H2_mother'].fillna(0)
    # replace the nan values in the motif_ids columns with empty string
    merged_df['motif_ids_H1_father'] = merged_df['motif_ids_H1_father'].fillna('')
    merged_df['motif_ids_H2_father'] = merged_df['motif_ids_H2_father'].fillna('')
    merged_df['motif_ids_H1_mother'] = merged_df['motif_ids_H1_mother'].fillna('')
    merged_df['motif_ids_H2_mother'] = merged_df['motif_ids_H2_mother'].fillna('')
    merged_df['motif_ids_H1_child'] = merged_df['motif_ids_H1_child'].fillna('')
    merged_df['motif_ids_H2_child'] = merged_df['motif_ids_H2_child'].fillna('')
    merged_df['CN_ref_father'] = merged_df['CN_ref_father'].fillna(0)
    merged_df['CN_ref_mother'] = merged_df['CN_ref_mother'].fillna(0)
    merged_df['CN_ref_child'] = merged_df['CN_ref_child'].fillna(0)
    merged_df['CN_H1_child'] = merged_df['CN_H1_child'].fillna(0)
    merged_df['CN_H2_child'] = merged_df['CN_H2_child'].fillna(0)

    merged_df['sex_child'] = child_sex

    merged_df['mendelian_consistent_0'] = merged_df.apply(lambda row: is_mendelian_consistent_based_on_motif(row), axis=1)
    merged_df['mendelian_consistent_1'] = merged_df.apply(lambda row: is_mendelian_consistent_based_on_motif(row), axis=1)
    
    #merged_df['mendelian_consistent_motif'] = merged_df.apply(lambda row: is_mendelian_consistent_based_on_motif(row, tolerance=0), axis=1)
    merged_df['tool'] = tool
    merged_df['max_motif_size'] = merged_df['motifs_child'].apply(max_motif_size)
    merged_df['region'] = merged_df['chrom'] + ':' + merged_df['start'].astype(str) + '-' + merged_df['end'].astype(str)

    return merged_df


def max_motif_size(motifs):
    if pd.isna(motifs):
        return 0  
    return max(len(i) for i in motifs)


def print_mendalian_consistency_per_motif_size(df):
    for i in sorted(df['max_motif_size'].unique()):
        print("Mendelian consistency for motif size ", i, " is " , df[df['max_motif_size'] == i]['mendelian_consistent'].sum()/df[df['max_motif_size'] == i].shape[0])
        print(df[df['max_motif_size'] == i]['mendelian_consistent'].value_counts())




def check_agreement(row):
    # Extract haplotype copy numbers for the three tools (tandem_twister, TRGT, vamos)
    twister = abs(row['CN_H1_tandem_twister']+row['CN_H2_tandem_twister'])
    trgt = abs(row['CN_H1_TRGT']+row['CN_H2_TRGT'])
    vamos = abs(row['CN_H1_vamos']+row['CN_H2_vamos'])
    
    # Check agreement between the tools considering haplotype swaps
    if twister == trgt == vamos:
        return 'all_agree'
    elif twister == trgt and twister != vamos:
        return 'vamos_disagree'
    elif twister == vamos and twister != trgt:
        return 'TRGT_disagree'
    elif trgt == vamos and trgt != twister:
        return 'tandem_twister_disagree'
    else:
        return 'all_disagree'

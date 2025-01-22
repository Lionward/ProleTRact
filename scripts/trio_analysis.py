import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from parse_vcf_files import parse_vcf_tandemtwister

from parse_vcf_files import parse_asm_vcf_tandemtwister

# Example DataFrame
# data = {
#     'chrom': ['chr1', 'chr1', 'chr2', 'chr2', 'chr3', 'chr3'],
#     'start': [1000, 2000, 1500, 2500, 3000, 4000],
#     'inheritance': [1, -1, 3, -3, 2, -2],
#     'child_hap': ['h1', 'h2', 'h1', 'h2', 'h1', 'h2']
# }
# trio_df = pd.DataFrame(data)

# def compute_trio_inheritance_df(child_h1_df, child_h2_df, mother_h1_df, mother_h2_df, father_h1_df, father_h2_df):
def compute_trio_inheritance_df(trio_dict):
    '''
    input: trio_dict
    a dictionary including these vcf files:
    child_h1_file, child_h2_file, mother_h1_file, mother_h2_file, father_h1_file, father_h2_file

    parsing files:
    the parsed files are stored in a dict with these dataframes:
    child_h1_df, child_h2_df, mother_h1_df, mother_h2_df, father_h1_df, father_h2_df


    '''
    
    # parsing data frames and merging them into trio_df
    
    child_h1_df = parse_asm_vcf_tandemtwister(pysam.VariantFile(trio_dict['child_h1_file']))[['chrom', 'start', 'end', 'motifs', 'CN']]
    child_h1_df.rename(columns={'CN': 'CN_h1'}, inplace=True)
    trio_df = child_h1_df

    child_h2_df = parse_asm_vcf_tandemtwister(pysam.VariantFile(trio_dict['child_h2_file']))[['chrom', 'start', 'end', 'CN']]
    child_h2_df.rename(columns={'CN': 'CN_h2'}, inplace=True)
    trio_df = pd.merge(trio_df, child_h2_df, on=['chrom', 'start', 'end'], how='inner')


    mother_h1_df = parse_asm_vcf_tandemtwister(pysam.VariantFile(trio_dict['mother_h1_file']))[['chrom', 'start', 'end', 'CN']]
    mother_h1_df.rename(columns={'CN': 'CN_mother_h1'}, inplace=True)
    trio_df = pd.merge(trio_df, mother_h1_df, on=['chrom', 'start', 'end'], how='inner')

    mother_h2_df = parse_asm_vcf_tandemtwister(pysam.VariantFile(trio_dict['mother_h2_file']))[['chrom', 'start', 'end', 'CN']]
    mother_h2_df.rename(columns={'CN': 'CN_mother_h2'}, inplace=True)
    trio_df = pd.merge(trio_df, mother_h2_df, on=['chrom', 'start', 'end'], how='inner')

    father_h1_df = parse_asm_vcf_tandemtwister(pysam.VariantFile(trio_dict['father_h1_file']))[['chrom', 'start', 'end', 'CN']]
    father_h1_df.rename(columns={'CN': 'CN_father_h1'}, inplace=True)
    trio_df = pd.merge(trio_df, father_h1_df, on=['chrom', 'start', 'end'], how='inner')

    father_h2_df = parse_asm_vcf_tandemtwister(pysam.VariantFile(trio_dict['father_h2_file']))[['chrom', 'start', 'end', 'CN']]
    father_h2_df.rename(columns={'CN': 'CN_father_h2'}, inplace=True)
    trio_df = pd.merge(trio_df, father_h2_df, on=['chrom', 'start', 'end'], how='inner')


    # analyzing trio_df, computing inheritance values

    # Filter out homozygous rows
    trio_df = trio_df[trio_df['CN_h1'] != trio_df['CN_h2']]
    trio_df = trio_df[trio_df['CN_mother_h1'] != trio_df['CN_mother_h2']]
    trio_df = trio_df[trio_df['CN_father_h1'] != trio_df['CN_father_h2']]

    # reduced to almost 7 % of the rows ....

    # assign parent to haplotype in each chromosome
    trio_df['d_h1_mother'] = trio_df.apply(lambda row: min(abs(row['CN_h1'] - row['CN_mother_h1']), abs(row['CN_h1'] - row['CN_mother_h2'])), axis=1)
    trio_df['d_h2_mother'] = trio_df.apply(lambda row: min(abs(row['CN_h2'] - row['CN_mother_h1']), abs(row['CN_h2'] - row['CN_mother_h2'])), axis=1)
    trio_df['d_h1_father'] = trio_df.apply(lambda row: min(abs(row['CN_h1'] - row['CN_father_h1']), abs(row['CN_h1'] - row['CN_father_h2'])), axis=1)
    trio_df['d_h2_father'] = trio_df.apply(lambda row: min(abs(row['CN_h2'] - row['CN_father_h1']), abs(row['CN_h2'] - row['CN_father_h2'])), axis=1)


    # Group by 'chrom' and sum the relevant columns
    sums = trio_df.groupby('chrom')[['d_h1_mother', 'd_h2_mother', 'd_h1_father', 'd_h2_father']].agg('sum')

    # Add the comparison to assign a parent based on the summed values
    sums['assign_parent_h1'] = sums.apply(
        lambda row: 'father' if (row['d_h1_mother'] + row['d_h2_father']) > (row['d_h1_father'] + row['d_h2_mother']) else 'mother', 
        axis=1
    )

    trio_df = trio_df.drop(columns=['d_h1_mother', 'd_h2_mother', 'd_h1_father', 'd_h2_father'])

    # # Merge the 'assign_parent_h1' back into the original DataFrame
    trio_df = trio_df.merge(sums[['assign_parent_h1']], on='chrom', how='left')

    def compute_inheritance(row, off_by=1):
        inheritance_h1, inheritance_h2 = None, None

        if row['assign_parent_h1'] == 'mother':
            # Compare h1 with mother and h2 with father
            if min(abs(row['CN_h1'] - row['CN_mother_h1']), abs(row['CN_h1'] - row['CN_mother_h2'])) <= off_by:
                inheritance_h1 = 'mother_h1' if abs(row['CN_h1'] - row['CN_mother_h1']) < abs(row['CN_h1'] - row['CN_mother_h2']) else 'mother_h2'
            else:
                inheritance_h1 = 'mother_denovo'

            if min(abs(row['CN_h2'] - row['CN_father_h1']), abs(row['CN_h2'] - row['CN_father_h2'])) <= off_by:
                inheritance_h2 = 'father_h1' if abs(row['CN_h2'] - row['CN_father_h1']) < abs(row['CN_h2'] - row['CN_father_h2']) else 'father_h2'
            else:
                inheritance_h2 = 'father_denovo'
        else:
            # Compare h2 with mother and h1 with father
            if min(abs(row['CN_h2'] - row['CN_mother_h1']), abs(row['CN_h2'] - row['CN_mother_h2'])) <= off_by:
                inheritance_h2 = 'mother_h1' if abs(row['CN_h2'] - row['CN_mother_h1']) < abs(row['CN_h2'] - row['CN_mother_h2']) else 'mother_h2'
            else:
                inheritance_h2 = 'mother_denovo'

            if min(abs(row['CN_h1'] - row['CN_father_h1']), abs(row['CN_h1'] - row['CN_father_h2'])) <= off_by:
                inheritance_h1 = 'father_h1' if abs(row['CN_h1'] - row['CN_father_h1']) < abs(row['CN_h1'] - row['CN_father_h2']) else 'father_h2'
            else:
                inheritance_h1 = 'father_denovo'

        return pd.Series({'inheritance_h1': inheritance_h1, 'inheritance_h2': inheritance_h2})

    # Apply the function row-wise and assign results to new columns
    trio_df[['inheritance_h1', 'inheritance_h2']] = trio_df.apply(compute_inheritance, axis=1)

    trio_df = trio_df[trio_df['chrom'] != 'chrX']

    # Custom chromosome order
    chrom_order = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
        "chr18", "chr19", "chr20", "chr21", "chr22"
    ]

    # Convert 'chrom' column to a categorical type with the custom order
    trio_df['chrom'] = pd.Categorical(trio_df['chrom'], categories=chrom_order, ordered=True)

    return trio_df


def plot_trio_inheritance(trio_df, alpha=0.2):
    '''
    input trio_df has these columns:
    chrom	start	end	motifs	CN_h1	CN_h2	CN_mother_h1	CN_mother_h2	CN_father_h1	CN_father_h2	inheritance_h1	inheritance_h2
    inheritance_h1 and inheritance_h2 \in ["mother_denovo", "mother_h2", "mother_h1", "father_h1", "father_h2", "father_denovo"]
    '''

    # Melt the dataframe to reshape it for plotting
    plot_df = trio_df.melt(
        id_vars=["chrom", "start"],
        value_vars=["inheritance_h1", "inheritance_h2"],
        var_name="inheritance_type",  # This column will specify 'inheritance_h1' or 'inheritance_h2'
        value_name="inheritance"      # This column will hold the actual inheritance values
    )


    # Ensure correct order of categories
    inheritance_order = ["mother_denovo", "mother_h2", "mother_h1", "father_h1", "father_h2", "father_denovo"]
    plot_df['inheritance'] = pd.Categorical(plot_df['inheritance'], categories=inheritance_order, ordered=True)

    # Set Seaborn style
    sns.set(style="whitegrid")

    # Create the scatter plot using FacetGrid
    facet = sns.FacetGrid(plot_df, col="chrom", col_wrap=5, height=3, sharex=False, sharey=True)
    facet.map_dataframe(
        sns.scatterplot,
        x="start",
        y="inheritance",   # Now refers to the melted column
        hue="inheritance_type",  # Differentiates between 'inheritance_h1' and 'inheritance_h2'
        palette={"inheritance_h1": "red", "inheritance_h2": "blue"},
        edgecolor=None,
        alpha=alpha
    )

    # Adjust x-axis labels to vertical
    # facet.set_xticklabels(rotation=90)

    # Add legend and adjust layout
    facet.add_legend(title="Inheritance Type")
    facet.set_axis_labels("Tandem Repeat ID (M)", "Inherited Haplotype Allele")
    facet.set_titles(col_template="{col_name}")

    # Adjust legend to improve visibility
    # for legend_artist in facet._legend.legendHandles:
    #     legend_artist.set_alpha(1)  # Set alpha to 1 for legend colors

    # Format the x-axis to show values in k-scale
    for ax in facet.axes.flatten():
        if ax is not None:
            ax.set_xticklabels([f"{int(x)/1000000:.1f}M" for x in ax.get_xticks()])

    # Improve spacing between subplots
    plt.subplots_adjust(top=0.9, hspace=0.3)
    
    return plt
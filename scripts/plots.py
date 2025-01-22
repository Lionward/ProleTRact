import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

def plot_heatmap_with_histograms(df, seq_sim_col='seq_sim', len_sim_col='len_sim', bins=30):
    """
    Plots a heatmap scatter plot with histograms on the sides for seq_sim and len_sim columns from a dataframe.

    Parameters:
        df (pd.DataFrame): DataFrame containing the data.
        seq_sim_col (str): Column name for x-axis values (default is 'seq_sim').
        len_sim_col (str): Column name for y-axis values (default is 'len_sim').
        bins (int): Number of bins for the histograms (default is 30).
    """
    
    # Set up the figure with subplots
    fig = plt.figure(figsize=(10, 10))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)

    # Main scatter heatmap
    ax_main = fig.add_subplot(grid[1:4, 0:3])
    scatter = ax_main.hexbin(df[seq_sim_col], df[len_sim_col], gridsize=50, cmap='coolwarm')
    ax_main.set_xlabel(seq_sim_col)
    ax_main.set_ylabel(len_sim_col)

    # Create color bar for the heatmap
    cbar = plt.colorbar(scatter, ax=ax_main, fraction=0.036, pad=0.04)
    cbar.set_label('Counts')

    # Create colormap for the histograms
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", "yellow", "red"])

    # Histogram for the x-axis
    ax_hist_x = fig.add_subplot(grid[0, 0:3], sharex=ax_main)
    counts, edges, _ = ax_hist_x.hist(df[seq_sim_col], bins=bins, color='grey', edgecolor='black')
    ax_hist_x.cla()  # Clear the axis to prepare for the colored histogram
    sns.histplot(df[seq_sim_col], bins=bins, color='gray', edgecolor='black', ax=ax_hist_x, kde=False)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    ax_hist_x.bar(bin_centers, counts, width=np.diff(edges), color=cmap((counts - counts.min()) / (counts.max() - counts.min())))
    ax_hist_x.set_ylabel('Frequency')
    ax_hist_x.set_title('Histogram of ' + seq_sim_col)
    ax_hist_x.set_yticks([])

    # Vertical histogram for the y-axis
    ax_hist_y = fig.add_subplot(grid[1:4, 3], sharey=ax_main)
    counts, edges, _ = ax_hist_y.hist(df[len_sim_col], bins=bins, orientation='horizontal', color='grey', edgecolor='black')
    ax_hist_y.cla()  # Clear the axis to prepare for the colored histogram
    sns.histplot(df[len_sim_col], bins=bins, color='gray', edgecolor='black', ax=ax_hist_y, kde=False)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    ax_hist_y.barh(bin_centers, counts, height=np.diff(edges), color=cmap((counts - counts.min()) / (counts.max() - counts.min())))
    ax_hist_y.set_xlabel('Frequency')
    ax_hist_y.set_title('Histogram of ' + len_sim_col, pad=20)
    ax_hist_y.set_xticks([])
    ax_hist_y.invert_xaxis()  # To align it properly with the y axis on the left.

    # Optional: Adjust the limits to make sure histograms don't overlap
    ax_main.set_xlim(df[seq_sim_col].min(), df[seq_sim_col].max())
    ax_main.set_ylim(df[len_sim_col].min(), df[len_sim_col].max())

    return plt

# Example usage:
# df = pd.DataFrame({
#     'seq_sim': np.random.rand(1000),
#     'len_sim': np.random.rand(1000)
# })

# print(df)

# heatmap_with_histograms(df).show()



def plot_karyotype_trio_hap_blocks(trio_df, output=None, title="Haplotype blocks"):

    # Example DataFrame with similar structure
    # data = {
    #     "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr4", "chr5"],
    #     "start": [5000000, 15000000, 3000000, 12000000, 1000000, 4000000, 2000000],
    #     "inheritance_h1": ["mother_h1", "father_h1", "mother_h1", "mother_h1", "father_h1", "mother_h1", "father_h1"],
    #     "inheritance_h2": ["mother_h2", "father_h2", "mother_h2", "father_h2", "mother_h2", "mother_h2", "father_h2"]
    # }

    # ######## trio_df = pd.DataFrame(data)

    # Sort chromosomes
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    trio_df["chrom"] = pd.Categorical(trio_df["chrom"], categories=chrom_order, ordered=True)
    trio_df = trio_df.sort_values(by=["chrom", "start"])

    # Haplotype color mapping
    haplotype_color_map = {"mother_h1": "red", "father_h1": "red", "mother_h2": "blue", "father_h2": "blue"}

    # Chromosome lengths for illustrative purposes
    chromosome_lengths = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468
    }

    # Reduce matplotlib calls
    fig, ax = plt.subplots(figsize=(12, 8))

    # Compute bar positions
    chrom_y_positions = {chrom: idx * 3 for idx, chrom in enumerate(chrom_order[::-1])}  # Reverse order for top-to-bottom

    # Draw karyotype bars
    for chrom, y_pos in chrom_y_positions.items():
        if chrom in chromosome_lengths:
            length = chromosome_lengths[chrom]
            ax.hlines(y=y_pos + 0.6, xmin=0, xmax=length, color="black", linewidth=5, alpha=0.5)
            ax.hlines(y=y_pos - 0.6, xmin=0, xmax=length, color="black", linewidth=5, alpha=0.5)

    # Create subsets for each chrom to speed up line plotting
    for chrom, group in trio_df.groupby("chrom"):
        if chrom not in chrom_y_positions:
            continue
        y_pos = chrom_y_positions[chrom]
        h1_data = group[group["inheritance_h1"].isin(haplotype_color_map.keys())]
        h2_data = group[group["inheritance_h2"].isin(haplotype_color_map.keys())]

        # Plot all inheritance_h1 in one call
        ax.vlines(
            x=h1_data["start"],
            ymin=y_pos + 0.6 - 0.2,
            ymax=y_pos + 0.6 + 0.2,
            color=h1_data["inheritance_h1"].map(haplotype_color_map),
            linewidth=2
        )

        # Plot all inheritance_h2 in one call
        ax.vlines(
            x=h2_data["start"],
            ymin=y_pos - 0.6 - 0.2,
            ymax=y_pos - 0.6 + 0.2,
            color=h2_data["inheritance_h2"].map(haplotype_color_map),
            linewidth=2
        )

        # Add haplotype labels (h1 and h2 on the left of bars)
        ax.text(-3e6, y_pos + 0.6, "h1", ha="right", va="center", fontsize=7)
        ax.text(-3e6, y_pos - 0.6, "h2", ha="right", va="center", fontsize=7)


    # Format plot
    ax.set_yticks([y for y in chrom_y_positions.values()])
    ax.set_yticklabels(chrom_order[::-1], fontsize=8)
    ax.set_xlabel("Genomic Position (bp)", fontsize=10)
    ax.set_ylabel("Chromosomes", fontsize=10)
    ax.set_title(title, fontsize=12)

    # Add legend
    legend_patches = [
        mpatches.Patch(color="red", label="parent's H1"),
        mpatches.Patch(color="blue", label="parent's H2")
    ]

    ax.legend(handles=legend_patches, loc="lower center")

    # Save plot to file
    plt.tight_layout()

    if output != None:
        plt.savefig(output, dpi=300)
    # plt.show()

    return plt
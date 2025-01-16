import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors

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

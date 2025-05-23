import os
import re
from typing import Dict, List

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


def load_mutation_data(file_path: str, sheet_name: str) -> pd.DataFrame:
    """Loads the mutation data from an Excel file and calculates mutation-related columns."""
    df = pd.read_excel(file_path, sheet_name=sheet_name)

    # Calculate the ancestor_MUT_lineage_MUT_ratio
    df["ancestor_MUT_lineage_MUT_ratio"] = round(df["ancestor_MUT_reads"] / df["lineage_MUT_READS"], 2)

    # Define mutation type
    df["mutation_type"] = df["ancestor_MUT_lineage_MUT_ratio"].apply(lambda x: "de_novo" if x == 0 else "ancestor")

    return df


def extract_lineage_info(df: pd.DataFrame, ancestor_phage: str) -> pd.DataFrame:
    """Extract ancestor phage, phage lineage, and host bacteria from combined_id."""

    df['Ancestor Phage'] = df['combined_id'].apply(lambda x: x.split('.')[0])
    df['Phage Lineage'] = df['combined_id'].apply(lambda x: '.'.join(x.split('.')[1:])[:-4])

    # Regex pattern to extract suffix (e.g., P1L1, P2L5, etc.)
    suffix_regexp = r'P\d+L\d+$'

    def extract_host_bacteria(lineage):
        return re.sub(suffix_regexp, '', lineage)  # Removes only the exact suffix

    # Apply the function to extract the correct Host Bacteria name
    df['Host Bacteria'] = df['Phage Lineage'].apply(extract_host_bacteria)

    # Print first 20 rows to verify the results
    # print(df[['combined_id', 'Phage Lineage', 'Host Bacteria']].head(20))

    return df


def load_infectivity_data(file_path: str) -> pd.DataFrame:
    """Loads infectivity data, replaces negative values with 0, and creates phage_name."""
    index_df = pd.read_csv(file_path)

    # Replace negative values in "Index" with 0
    index_df["Index"] = index_df["Index"].apply(lambda x: max(0, x))

    # Create "phage_name" column
    index_df["phage_name"] = index_df["strain"] + index_df["phage"] + "L" + index_df["lineage"].astype(str)

    # print("\n🔹 Infectivity Data Before Merging:")
    # print(index_df.head())

    return index_df


def merge_infectivity_data(df: pd.DataFrame, index_df: pd.DataFrame) -> pd.DataFrame:
    """Merges the infectivity Index into the main df based on Phage Lineage matching phage_name."""
    # print("\n🔹 Mutation Data Before Merging:")
    # print(df.head())

    merged_df = df.merge(index_df[["phage_name", "Index"]], left_on="Phage Lineage", right_on="phage_name", how="left")

    # print("\n🔹 Merged Data After Merging Infectivity Index:")
    # print(merged_df.head())

    return merged_df


def parse_genbank(genbank_file: str) -> List[Dict]:
    """Parses a GenBank file to extract gene annotations."""
    genes = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start)
                end = int(feature.location.end)
                gene_name = feature.qualifiers.get("product", ["unknown"])[0]
                genes.append({"start": start, "end": end, "gene": gene_name})
    return genes


def plot_gene_map(ax, genes, gene_y):
    print("Plotting gene map at y=", gene_y)  # Debugging
    y_offset = 0.6

    arrow_body_width = 1.5  # thicker shaft
    arrow_head_length = 240  # same head
    arrow_head_width = 4  # same width as before, or tweak if needed

    for i, gene in enumerate(genes):

        ax.arrow(
            gene["start"], gene_y,
            gene["end"] - gene["start"], 0,
            width=arrow_body_width,
            head_length=arrow_head_length,
            head_width=arrow_head_width,
            length_includes_head=True,
            fc='gray', ec='black', zorder=50
        )

        ax.text((gene["start"] + gene["end"]) / 2, gene_y + y_offset * ((-1) ** i),
                gene["gene"], ha='center', fontsize=28, rotation=90, va='bottom', zorder=50)  # ✅ Ensure text is on top


def plot_ancestor_line(ax, ancestor_y, unique_reference_positions, ancestor_phage, mutation_colors):
    print(f"Plotting ancestor line at y={ancestor_y}")

    ax.plot([0, 6034], [ancestor_y, ancestor_y], linestyle='-', color='black', alpha=0.8, linewidth=3, zorder=50)

    ref_colors = unique_reference_positions['REF'].map(mutation_colors).fillna('black')

    ax.scatter(unique_reference_positions['POS'], [ancestor_y] * len(unique_reference_positions),
               c=ref_colors, edgecolors='black', alpha=0.8, s=300, linewidths=0.3, zorder=50)

    ax.text(-600, ancestor_y, ancestor_phage, va='center', fontsize=28, fontweight='bold', ha='right', zorder=50)


def plot_phage_mutations(ax, df, lineage_map, first_lineage_per_host):
    """Plots mutation lines and points for each phage lineage using the mutation_type column for color."""
    mutation_colors = {"de_novo": "red", "ancestor": "black"}

    print("Unique Host Names:", df['Host Bacteria'].unique())

    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        ax.plot([0, 6034], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)

        # Assign colors based on mutation_type
        colors = data["mutation_type"].map(mutation_colors)

        y_pos = lineage_map[lineage]  # Make sure it uses the new spaced y-positions

        for pos, color in zip(data['POS'], colors):
            ax.vlines(x=pos, ymin=y_pos - 0.4, ymax=y_pos + 0.4, color=color, linewidth=2, alpha=0.9)
        if lineage in first_lineage_per_host:
            ax.text(-600, y_pos, first_lineage_per_host[lineage], va='center', fontsize=28, ha='right')

    # Ensure consistent font size for x and y axis ticks
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)


def plot_stacked_mutation_histogram(ax, lineage_map, mutation_counts, de_novo_counts):
    y_positions = [lineage_map[l] for l in lineage_map]
    total_values = [mutation_counts.get(l, 0) for l in lineage_map]
    de_novo_values = [de_novo_counts.get(l, 0) for l in lineage_map]
    ancestor_values = [t - d for t, d in zip(total_values, de_novo_values)]

    ax.barh(y_positions, de_novo_values, color='red', alpha=0.8, height=0.5,
            edgecolor='black', label='De Novo')
    ax.barh(y_positions, ancestor_values, left=de_novo_values, color='gray',
            alpha=0.6, height=0.5, edgecolor='black', label='Ancestor')

    # Add labels
    for y, t in zip(y_positions, total_values):
        ax.text(t + 0.5, y, str(t), va='center', fontsize=20)

    ax.set_xlabel("Mutation Count", fontsize=28, fontweight='bold')
    ax.set_ylabel("Phage Lineage", fontsize=28, fontweight='bold')
    ax.set_title("Total vs. De Novo Mutations", fontsize=24, fontweight='bold')

    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(max(1, max(total_values) // 5)))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    ax.legend(fontsize=20)


def plot_mutations(df: pd.DataFrame, genes: List[Dict], ancestor_phage: str, output_path: str, keep_lineage_name: bool = False) -> None:
    """Main function to plot the mutations along the genome with histograms and a heatmap."""

    lineages = list(df['Phage Lineage'].unique())  # Do NOT include the ancestor

    # Adjust y-axis spacing to align all elements properly
    host_bacteria_groups = df.groupby(df['Host Bacteria'].str.strip())["Phage Lineage"].unique()

    # Initialize lineage map with extra spacing
    lineage_map = {}

    y_position = 0
    extra_space = 2  # Extra spacing after each Host Bacteria group

    for host, lineages_in_host in host_bacteria_groups.items():
        for lineage in reversed(lineages_in_host):
            lineage_map[lineage] = y_position
            y_position += 1
        y_position += extra_space  # Extra spacing

    first_lineage_per_host = {}
    for host, lineages_in_host in host_bacteria_groups.items():
        if len(lineages_in_host) > 0:
            first_lineage_per_host[lineages_in_host[0]] = host

    gene_y = max(lineage_map.values()) + 3

    suffix_regexp = r'P\d+L\d+$'
    lineage_labels = {}

    for lineage in lineage_map:
        if keep_lineage_name:
            lineage_labels[lineage] = lineage
        else:
            lineage_labels[lineage] = re.sub(suffix_regexp, '', lineage)

    df['Lineage Order'] = df['Phage Lineage'].map(lineage_map)

    mutation_counts = df.groupby('Phage Lineage')['POS'].count().to_dict()

    # Mutation counts for de_novo only
    de_novo_counts = df[df["mutation_type"] == "de_novo"].groupby("Phage Lineage")["POS"].count().to_dict()

    fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [3, 1]},
                            figsize=(30, len(lineages) * 0.6 + 10))
    ax, ax_hist = axs

    ax_hist.sharey(ax)

    ax.set_ylim(-5, max(lineage_map.values()) + 3)

    # 🔹 **Fix subplot positions**
    fig.subplots_adjust(left=0.1, right=0.95, top=1.0, bottom=0.05, wspace=0.1)

    ax_hist.set_position([ax_hist.get_position().x0, ax.get_position().y0,
                          ax_hist.get_position().width, ax.get_position().height])

    # 🔹 Ensure Phage Mutations are Drawn
    plot_phage_mutations(ax, df, lineage_map, first_lineage_per_host)

    # 🔹 Ensure Histograms are Drawn Correctly
    plot_stacked_mutation_histogram(ax_hist, lineage_map, mutation_counts, de_novo_counts)

    # ✅ Plot the reference genome and gene map last to ensure they are visible
    plot_gene_map(ax, genes, gene_y)

    # ✅ Re-set y-axis to ensure visibility
    ax.set_ylim(-5, gene_y + 5)
    ax.set_yticks([])

    print("Gene Y:", gene_y)
    # print("Ancestor Y:", ancestor_y)
    print("Lineage Map:", lineage_map)

    # 🔹 **Ensure output directory exists and save**
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == "__main__":
    # File paths
    wd = "/mnt/c/crassvirales/wallapat_project"
    file_path = f"{wd}/20250213_REF_MUT_analysis.xlsx"
    infectivity_file = f"{wd}/infectivity_index_P1-P2.csv"

    results = f"{wd}/results"
    figures = f"{results}/figures"

    # Load infectivity index data
    index_df = load_infectivity_data(infectivity_file)

    phages = ('P1L1', 'P2L1')
    for ancestor_phage in phages:
        genbank_file = f"{wd}/{ancestor_phage}_phold_annotation.gbk"
        sheet_name = f"{ancestor_phage}_filter"

        output_file = f"{figures}/{ancestor_phage}_phage_mutations_with_genes_de_novo.png"

        # Load mutation data
        df = load_mutation_data(file_path, sheet_name)
        df = extract_lineage_info(df, ancestor_phage=ancestor_phage)

        # Merge infectivity index
        df = merge_infectivity_data(df, index_df)

        # Load gene annotations
        genes = parse_genbank(genbank_file)

        # Plot mutations
        plot_mutations(df, genes, f'{ancestor_phage}_reference', output_file)

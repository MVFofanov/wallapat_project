import os
from typing import Dict, List

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
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
    df['Host Bacteria'] = df['Phage Lineage'].apply(lambda x: x.rsplit(ancestor_phage, 1)[0])
    return df


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
    """Plots the gene annotation map."""
    y_offset = 0.6
    for i, gene in enumerate(genes):
        ax.arrow(gene["start"], gene_y, gene["end"] - gene["start"], 0, head_width=0.5, head_length=100,
                 fc='gray', ec='black', length_includes_head=True)
        ax.text((gene["start"] + gene["end"]) / 2, gene_y + y_offset * ((-1) ** i),
                gene["gene"], ha='center', fontsize=24, rotation=90, va='bottom')


def plot_ancestor_line(ax, ancestor_y, unique_reference_positions, ancestor_phage, mutation_colors):
    """Plots the reference genome line with reference nucleotide positions."""
    ax.plot([0, 6034], [ancestor_y, ancestor_y], linestyle='-', color='black', alpha=0.6, linewidth=1.5)

    # Ensure REF values are mapped properly, replace NaN with default color
    ref_colors = unique_reference_positions['REF'].map(mutation_colors).fillna('gray')

    ax.scatter(unique_reference_positions['POS'], [ancestor_y] * len(unique_reference_positions),
               c=ref_colors, edgecolors='black', alpha=0.8, s=200, linewidths=0.2)

    ax.text(-600, ancestor_y, ancestor_phage, va='center', fontsize=24, fontweight='bold', ha='right')


def plot_phage_mutations(ax, df, lineage_map):
    """Plots mutation lines and points for each phage lineage using the mutation_type column for color."""
    mutation_colors = {"de_novo": "red", "ancestor": "black"}

    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        ax.plot([0, 6034], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)

        # Assign colors based on mutation_type
        colors = data["mutation_type"].map(mutation_colors)

        ax.scatter(data['POS'], [y_pos] * len(data), c=colors, label=lineage, alpha=0.8,
                   edgecolors='black', s=200, linewidths=0.2)
        ax.text(-600, y_pos, lineage, va='center', fontsize=24, ha='right')

    # Ensure consistent font size for x and y axis ticks
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)


def plot_mutation_histogram(ax_hist, lineage_map, mutation_counts):
    """Plots the histogram of mutation counts per lineage."""
    y_positions = [lineage_map[l] for l in lineage_map.keys()]
    hist_values = [mutation_counts[l] for l in lineage_map.keys()]

    ax_hist.barh(y_positions, hist_values, color='gray', alpha=0.6, height=0.4, align='center', edgecolor='black')

    for y, count in zip(y_positions, hist_values):
        ax_hist.text(count + 1, y, str(count), va='center', fontsize=14, fontweight='bold')

    ax_hist.set_xlabel("Total Mutations", fontsize=28, fontweight='bold')
    ax_hist.set_ylabel("Phage Lineage", fontsize=28, fontweight='bold')
    ax_hist.set_title("Mutation Count per Lineage", fontsize=24, fontweight='bold')

    ax_hist.xaxis.set_tick_params(labelsize=24, labelbottom=True)
    ax_hist.yaxis.set_tick_params(labelsize=24)

    ax_hist.xaxis.set_major_locator(ticker.MultipleLocator(max(1, max(hist_values) // 5)))
    ax_hist.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax_hist.grid(axis="x", linestyle="--", alpha=0.5)


def plot_de_novo_histogram(ax_de_novo, lineage_map, de_novo_counts):
    """Plots a histogram of de_novo mutations per lineage."""
    y_positions = [lineage_map[l] for l in lineage_map.keys()]
    hist_values = [de_novo_counts.get(l, 0) for l in lineage_map.keys()]  # Default to 0 if no de_novo mutations

    ax_de_novo.barh(y_positions, hist_values, color='red', alpha=0.6, height=0.4, align='center', edgecolor='black')

    for y, count in zip(y_positions, hist_values):
        ax_de_novo.text(count + 1, y, str(count), va='center', fontsize=14, fontweight='bold')

    ax_de_novo.set_xlabel("De Novo Mutations", fontsize=28, fontweight='bold')
    ax_de_novo.set_title("De Novo Mutation Count", fontsize=24, fontweight='bold')

    ax_de_novo.xaxis.set_tick_params(labelsize=24, labelbottom=True)
    ax_de_novo.yaxis.set_tick_params(labelsize=24)

    ax_de_novo.xaxis.set_major_locator(ticker.MultipleLocator(max(1, max(hist_values) // 5)))
    ax_de_novo.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax_de_novo.grid(axis="x", linestyle="--", alpha=0.5)

    # Remove y-axis labels on this histogram (they're on the left plot)
    ax_de_novo.set_yticks([])
    ax_de_novo.set_yticklabels([])


def plot_mutations(df: pd.DataFrame, genes: List[Dict], ancestor_phage: str, output_path: str) -> None:
    """Main function to plot the mutations along the genome with two histograms (all & de novo mutations)."""

    unique_reference_positions = df[['POS', 'REF']].drop_duplicates()
    lineages = list(df['Phage Lineage'].unique())
    lineages.insert(0, ancestor_phage)
    lineage_map = {lineage: i for i, lineage in enumerate(reversed(lineages))}
    df['Lineage Order'] = df['Phage Lineage'].map(lineage_map)

    mutation_counts = df.groupby('Phage Lineage')['POS'].count().to_dict()
    mutation_counts[ancestor_phage] = len(unique_reference_positions)

    # Mutation counts for de_novo only
    de_novo_counts = df[df["mutation_type"] == "de_novo"].groupby("Phage Lineage")["POS"].count().to_dict()

    fig, axs = plt.subplots(
        nrows=1, ncols=3, gridspec_kw={'width_ratios': [3, 1, 1]}, figsize=(30, len(lineages) * 0.6 + 6)
    )

    ax, ax_hist, ax_de_novo = axs
    ax_hist.sharey(ax)  # Align first histogram with left plot
    ax_de_novo.sharey(ax)  # Align second histogram with left plot

    plot_gene_map(ax, genes, len(lineages) + 1)
    plot_ancestor_line(ax, lineage_map[ancestor_phage], unique_reference_positions, ancestor_phage, {'A': 'gray'})
    plot_phage_mutations(ax, df, lineage_map)
    plot_mutation_histogram(ax_hist, lineage_map, mutation_counts)
    plot_de_novo_histogram(ax_de_novo, lineage_map, de_novo_counts)

    plt.subplots_adjust(wspace=0.2)  # Adjust spacing between plots

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=900)
    plt.close()


if __name__ == "__main__":
    # File path and sheet name (replace with actual path)
    wd = "/mnt/c/crassvirales/wallapat_project"
    file_path = f"{wd}/20250213_REF_MUT_analysis.xlsx"

    results = f"{wd}/results"
    figures = f"{results}/figures"

    phages = ('P1L1', 'P2L1')
    for ancestor_phage in phages:
        genbank_file = f"{wd}/{ancestor_phage}_phold_annotation.gbk"
        sheet_name = f"{ancestor_phage}_filter"

        output_file = f"{figures}/{ancestor_phage}_phage_mutations_with_genes_de_novo.png"

        df = load_mutation_data(file_path, sheet_name)
        df = extract_lineage_info(df, ancestor_phage=ancestor_phage)
        genes = parse_genbank(genbank_file)
        plot_mutations(df, genes, f'{ancestor_phage}_reference', output_file)

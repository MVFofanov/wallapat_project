import os
from typing import Dict, List

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


def load_mutation_data(file_path: str, sheet_name: str) -> pd.DataFrame:
    """Loads the mutation data from an Excel file."""
    return pd.read_excel(file_path, sheet_name=sheet_name)


def extract_lineage_info(df: pd.DataFrame, ancestor_phage: str) -> pd.DataFrame:
    """Extract ancestor phage, phage lineage, and host bacteria from combined_id."""
    df['Ancestor Phage'] = df['combined_id'].apply(lambda x: x.split('.')[0])
    df['Phage Lineage'] = df['combined_id'].apply(lambda x: '.'.join(x.split('.')[1:]))
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
                gene["gene"], ha='center', fontsize=8, rotation=45, va='bottom')


def plot_ancestor_line(ax, ancestor_y, unique_reference_positions, ancestor_phage, mutation_colors):
    """Plots the reference genome line with reference nucleotide positions."""
    ax.plot([0, 6034], [ancestor_y, ancestor_y], linestyle='-', color='black', alpha=0.6, linewidth=1.5)
    ax.scatter(unique_reference_positions['POS'], [ancestor_y] * len(unique_reference_positions),
               c=unique_reference_positions['REF'].map(mutation_colors), edgecolors='black', alpha=0.8, s=60,
               linewidths=0.2)
    ax.text(-600, ancestor_y, ancestor_phage, va='center', fontsize=10, fontweight='bold', ha='right')


def plot_phage_mutations(ax, df, lineage_map, mutation_colors):
    """Plots mutation lines and points for each phage lineage."""
    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        ax.plot([0, 6034], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)
        ax.scatter(data['POS'], [y_pos] * len(data), c=data['MUT'].map(mutation_colors), label=lineage, alpha=0.8,
                   edgecolors='black', s=60, linewidths=0.2)
        ax.text(-600, y_pos, lineage, va='center', fontsize=10, fontweight='bold', ha='right')


def plot_mutation_histogram(ax_hist, lineage_map, mutation_counts):
    """Plots the histogram of mutation counts per lineage."""
    y_positions = [lineage_map[l] for l in lineage_map.keys()]
    hist_values = [mutation_counts[l] for l in lineage_map.keys()]
    ax_hist.barh(y_positions, hist_values, color='gray', alpha=0.6, height=0.4, align='center', edgecolor='black')
    for y, count in zip(y_positions, hist_values):
        ax_hist.text(count + 1, y, str(count), va='center', fontsize=10)


def plot_stacked_barplot(ax_bar, df, mutation_colors):
    """Plots a stacked barplot for mutation positions, excluding the reference genome."""
    mutation_counts_per_pos = df.groupby(['POS', 'MUT']).size().unstack(fill_value=0)
    bottom_values = None
    for nucleotide in ['A', 'C', 'G', 'T']:
        if nucleotide in mutation_counts_per_pos.columns:
            ax_bar.bar(mutation_counts_per_pos.index, mutation_counts_per_pos[nucleotide], bottom=bottom_values,
                       color=mutation_colors[nucleotide], label=nucleotide)
            bottom_values = mutation_counts_per_pos[nucleotide] if bottom_values is None else bottom_values + mutation_counts_per_pos[nucleotide]
    ax_bar.set_xlabel("Genomic Position")
    ax_bar.set_ylabel("Mutation Frequency")
    ax_bar.set_title("Stacked Nucleotide Mutations per Position")
    ax_bar.legend()


def plot_mutations(df: pd.DataFrame, genes: List[Dict], ancestor_phage: str, output_path: str) -> None:
    """Main function to plot the mutations along the genome."""
    mutation_colors = {'A': '#377eb8', 'C': '#ff7f00', 'G': '#984ea3', 'T': '#a65628'}

    unique_reference_positions = df[['POS', 'REF']].drop_duplicates()
    lineages = list(df['Phage Lineage'].unique())
    lineages.insert(0, ancestor_phage)
    lineage_map = {lineage: i for i, lineage in enumerate(reversed(lineages))}
    df['Lineage Order'] = df['Phage Lineage'].map(lineage_map)

    mutation_counts = df.groupby('Phage Lineage')['POS'].count().to_dict()
    mutation_counts[ancestor_phage] = len(unique_reference_positions)

    fig, axs = plt.subplots(nrows=2, ncols=2, gridspec_kw={'height_ratios': [3, 1], 'width_ratios': [3, 1]},
                            figsize=(24, len(lineages) * 0.6 + 6), sharey='row')
    ax, ax_hist = axs[0]
    ax_bar, ax_empty = axs[1]

    gene_y = len(lineages) + 1
    plot_gene_map(ax, genes, gene_y)
    ancestor_y = lineage_map[ancestor_phage]
    plot_ancestor_line(ax, ancestor_y, unique_reference_positions, ancestor_phage, mutation_colors)
    plot_phage_mutations(ax, df, lineage_map, mutation_colors)
    plot_mutation_histogram(ax_hist, lineage_map, mutation_counts)
    plot_stacked_barplot(ax_bar, df[df['Phage Lineage'] != ancestor_phage], mutation_colors)

    ax.set_xlabel("Genomic Position")
    ax.set_yticks([])
    ax.set_xlim(-600, 6500)
    ax.set_title("Phage Lineage Mutations and Genomic Map")
    ax_hist.set_xlabel("Number of Mutations")
    ax_hist.set_ylabel("Phage Lineage")
    ax_hist.set_title("Mutation Count per Lineage")
    ax_bar.set_xlabel("Genomic Position")
    ax_bar.set_ylabel("Mutation Frequency")
    ax_bar.set_title("Stacked Nucleotide Mutations per Position")

    ax_empty.set_xticks([])
    ax_empty.set_yticks([])
    ax_empty.axis("off")

    ax.legend(
        handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=n) for n, c in
                 mutation_colors.items()], title="Mutation Type", loc='upper left', bbox_to_anchor=(1, 1))

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=900)
    plt.close()


if __name__ == "__main__":
    # File path and sheet name (replace with actual path)
    wd = "/mnt/c/crassvirales/wallapat_project"
    file_path = f"{wd}/20250213_REF_MUT_analysis.xlsx"
    ancestor_phage = 'P1L1'
    genbank_file = f"{wd}/{ancestor_phage}_phold_annotation.gbk"
    sheet_name = f"{ancestor_phage}_filter"

    results = f"{wd}/results"
    figures = f"{results}/figures"
    output_file = f"{figures}/phage_mutations_with_genes.png"

    df = load_mutation_data(file_path, sheet_name)
    df = extract_lineage_info(df, ancestor_phage=ancestor_phage)
    genes = parse_genbank(genbank_file)
    plot_mutations(df, genes, f'{ancestor_phage}_reference', output_file)

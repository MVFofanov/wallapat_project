import os
from typing import Dict, Tuple, List

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
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


def plot_mutations(df: pd.DataFrame, genes: List[Dict], ancestor_phage: str, output_path: str) -> None:
    """Plots the mutations along the genome for each phage lineage with a genomic map, reference genome line,
    and a horizontal histogram of mutation counts."""
    # Define colorblind-friendly mutation colors
    mutation_colors: Dict[str, str] = {
        'A': '#377eb8',  # Blue
        'C': '#ff7f00',  # Orange
        'G': '#984ea3',  # Purple
        'T': '#a65628'  # Brown
    }

    # Remove duplicate positions for reference genome
    unique_reference_positions = df[['POS', 'REF']].drop_duplicates()

    # Order lineages based on occurrence
    lineages = list(df['Phage Lineage'].unique())
    lineages.insert(0, ancestor_phage)  # Add ancestor phage as the first line
    lineage_map = {lineage: i for i, lineage in enumerate(reversed(lineages))}
    df['Lineage Order'] = df['Phage Lineage'].map(lineage_map)

    # Count mutations per lineage
    mutation_counts = df.groupby('Phage Lineage')['POS'].count().to_dict()
    mutation_counts[ancestor_phage] = len(unique_reference_positions)  # Updated count for ancestor

    fig, (ax, ax_hist) = plt.subplots(ncols=2, gridspec_kw={'width_ratios': [3, 1]},
                                      figsize=(24, len(lineages) * 0.6 + 3), sharey=True)  # Increased width

    # Plot gene map
    gene_y = len(lineages) + 1  # Position above mutations
    y_offset = 0.6  # Adjust vertical spacing for text labels

    for i, gene in enumerate(genes):
        ax.arrow(gene["start"], gene_y, gene["end"] - gene["start"], 0, head_width=0.5, head_length=100,
                 fc='gray', ec='black', length_includes_head=True)
        ax.text((gene["start"] + gene["end"]) / 2, gene_y + y_offset * ((-1) ** i),
                gene["gene"], ha='center', fontsize=8, rotation=45, va='bottom')

    # Plot ancestor genome line
    ancestor_y = lineage_map[ancestor_phage]
    ax.plot([0, 6034], [ancestor_y, ancestor_y], linestyle='-', color='black', alpha=0.6, linewidth=1.5)
    ax.scatter(unique_reference_positions['POS'], [ancestor_y] * len(unique_reference_positions),
               c=unique_reference_positions['REF'].map(mutation_colors), edgecolors='black', alpha=0.8, s=60,
               linewidths=0.2)
    ax.text(-600, ancestor_y, ancestor_phage, va='center', fontsize=10, fontweight='bold', ha='right')

    # Plot lines for each phage lineage
    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        ax.plot([0, 6034], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)
        ax.scatter(data['POS'], [y_pos] * len(data), c=data['MUT'].map(mutation_colors), label=lineage, alpha=0.8,
                   edgecolors='black', s=60, linewidths=0.2)
        ax.text(-600, y_pos, lineage, va='center', fontsize=10, fontweight='bold', ha='right')

    # Plot histogram of mutation counts
    y_positions = [lineage_map[l] for l in lineages]
    hist_values = [mutation_counts[l] for l in lineages]
    ax_hist.barh(y_positions, hist_values, color='gray', alpha=0.6, height=0.4, align='center')
    for y, count in zip(y_positions, hist_values):
        ax_hist.text(count + 1, y, str(count), va='center', fontsize=10)

    # Labels and aesthetics
    ax.set_xlabel("Genomic Position")
    ax.set_yticks(y_positions)
    ax.set_yticklabels([])
    ax.set_xlim(-600, 6500)  # moves start position for phage line name
    ax.set_title("Phage Lineage Mutations and Genomic Map")

    ax_hist.set_xlim(0, max(hist_values) * 1.1)
    ax_hist.set_xticks(range(0, max(hist_values) + 1, max(1, max(hist_values) // 5)))
    ax_hist.set_xlabel("Number of Mutations")
    ax_hist.set_yticks(y_positions)
    ax_hist.set_yticklabels([])
    ax_hist.set_title("Mutation Count")

    # Move legend outside of the plot
    ax.legend(
        handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=n) for n, c in
                 mutation_colors.items()], title="Mutation Type", loc='upper left', bbox_to_anchor=(1, 1))

    plt.tight_layout()

    # Ensure output directory exists
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

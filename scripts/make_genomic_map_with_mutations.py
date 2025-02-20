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

def plot_mutations(df: pd.DataFrame, genes: List[Dict], output_path: str) -> None:
    """Plots the mutations along the genome for each phage lineage with a genomic map."""
    # Define colorblind-friendly mutation colors
    mutation_colors: Dict[str, str] = {
        'A': '#377eb8',  # Blue
        'C': '#ff7f00',  # Orange
        'G': '#984ea3',  # Purple
        'T': '#a65628'  # Brown
    }

    # Order lineages based on occurrence
    lineages = df['Phage Lineage'].unique()
    lineage_map = {lineage: i for i, lineage in enumerate(reversed(lineages))}
    df['Lineage Order'] = df['Phage Lineage'].map(lineage_map)

    fig, ax = plt.subplots(figsize=(12, len(lineages) * 0.6 + 2))

    # Plot gene map
    gene_y = len(lineages) + 1  # Position above mutations
    for gene in genes:
        ax.arrow(gene["start"], gene_y, gene["end"] - gene["start"], 0, head_width=0.5, head_length=100,
                 fc='gray', ec='black', length_includes_head=True)
        ax.text((gene["start"] + gene["end"]) / 2, gene_y + 0.5, gene["gene"], ha='center', fontsize=8)

    # Plot lines for each phage lineage
    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        ax.plot([0, 6034], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)
        ax.scatter(data['POS'], [y_pos] * len(data), c=data['MUT'].map(mutation_colors), label=lineage,
                   edgecolors='black', s=60)
        ax.text(-600, y_pos, lineage, va='center', fontsize=10, fontweight='bold', ha='right')

    # Labels and aesthetics
    ax.set_xlabel("Genomic Position")
    ax.set_yticks([])
    # ax.set_xlim(-800, 6500)
    ax.set_xlim(-600, 6500) # moves start position for phage line name
    ax.set_title("Phage Lineage Mutations and Genomic Map")

    # Move legend outside of the plot
    ax.legend(
        handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=n) for n, c in
                 mutation_colors.items()], title="Mutation Type", loc='upper left', bbox_to_anchor=(1, 1))

    plt.tight_layout()

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight')
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
    plot_mutations(df, genes, output_file)

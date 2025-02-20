import os
from typing import Dict, Tuple

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

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


def plot_mutations(df: pd.DataFrame, output_path: str) -> None:
    """Plots the mutations along the genome for each phage lineage and saves the figure."""
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

    plt.figure(figsize=(12, len(lineages) * 0.6))

    # Plot lines for each phage lineage
    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        plt.plot([0, 6034], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)
        plt.scatter(data['POS'], [y_pos] * len(data), c=data['MUT'].map(mutation_colors), label=lineage,
                    edgecolors='black', s=60)
        plt.text(-600, y_pos, lineage, va='center', fontsize=10, fontweight='bold', ha='right', bbox=None)

    # Labels and aesthetics
    plt.xlabel("Genomic Position")
    plt.yticks([])
    plt.xlim(-800, 6500)
    plt.title("Phage Lineage Mutations")

    # Move legend outside of the plot
    plt.legend(
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
    sheet_name = f"{ancestor_phage}_filter"

    results = f"{wd}/results"
    figures = f"{results}/figures"
    output_file = f"{figures}/phage_mutations.png"

    df = load_mutation_data(file_path, sheet_name)
    df = extract_lineage_info(df, ancestor_phage=ancestor_phage)
    plot_mutations(df, output_file)

import os
import re
from typing import Dict, List, Union

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from Bio import SeqIO

matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


FUNCTION_COLORS = {
    "head and packaging": "#ff008d",
    "transcription regulation": "#ffe700",
    "lysis": "#001eff",
    "dna, rna and nucleotide metabolism": "#f000ff",
    "tail": "#74ee15",
    "connector": "#5A5A5A",
    "integration and excision": "#E0B0FF",
    "moron, auxiliary metabolic gene and host takeover": "#8900ff",
    "other": "#4deeea",
    "unknown": "#AAAAAA",
    "unknown function": "#AAAAAA"
}


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


def parse_genbank(genbank_file: str) -> List[Dict[str, str]]:
    """Parses a GenBank file to extract gene annotations and their functions."""
    genes = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start)
                end = int(feature.location.end)
                gene_name = feature.qualifiers.get("product", ["unknown"])[0]
                function = feature.qualifiers.get("function", ["unknown function"])[0]
                genes.append({
                    "start": start,
                    "end": end,
                    "gene": gene_name,
                    "function": function
                })
    return genes


def get_genome_length(genbank_file: str) -> int:
    """Returns the full genome length from a GenBank file."""
    for record in SeqIO.parse(genbank_file, "genbank"):
        return len(record.seq)


def plot_gene_map(ax: matplotlib.axes.Axes, genes: List[Dict[str, Union[str, int]]], gene_y: float) -> None:
    print("Plotting gene map at y=", gene_y)  # Debugging

    for i, gene in enumerate(genes):
        gene_start, gene_end = gene["start"], gene["end"]
        gene_mid = (gene_start + gene_end) / 2

        function = gene.get("function", "unknown function").strip().lower()

        # Normalize the function key to match dictionary
        color = FUNCTION_COLORS.get(function, FUNCTION_COLORS["unknown"])

        ax.arrow(
            gene_start, gene_y,
            gene_end - gene_start, 0,
            width=0.5,
            head_length=120,
            head_width=1,
            length_includes_head=True,
            fc=color, ec='black', zorder=50
        )

        ax.text(
            gene_mid, gene_y + 0.5,  # fixed offset
            gene["gene"], ha='center', fontsize=28, rotation=90,
            va='bottom', zorder=50
        )


def plot_phage_mutations(ax: matplotlib.axes.Axes,
                         df: pd.DataFrame,
                         lineage_map: Dict[str, int],
                         first_lineage_per_host: Dict[str, str],
                         genome_length: int
                         ) -> None:
    """Plots mutation lines and points for each phage lineage using the mutation_type column for color."""
    mutation_colors = {"de_novo": "red", "ancestor": "black"}

    print("Unique Host Names:", df['Host Bacteria'].unique())

    for lineage, data in df.groupby('Phage Lineage'):
        y_pos = lineage_map[lineage]
        ax.plot([0, genome_length], [y_pos, y_pos], linestyle='-', color='gray', alpha=0.5)

        # Assign colors based on mutation_type
        colors = data["mutation_type"].map(mutation_colors)

        y_pos = lineage_map[lineage]  # Make sure it uses the new spaced y-positions

        for pos, color in zip(data['POS'], colors):
            ax.vlines(x=pos, ymin=y_pos - 0.2, ymax=y_pos + 0.2, color=color, linewidth=2, alpha=0.9)
        if lineage in first_lineage_per_host:
            ax.text(-600, y_pos, first_lineage_per_host[lineage], va='center', fontsize=28, ha='right')

    # Ensure consistent font size for x and y axis ticks
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)


def plot_stacked_mutation_histogram(ax: matplotlib.axes.Axes,
                                    lineage_map: Dict[str, int],
                                    mutation_counts: Dict[str, int],
                                    de_novo_counts: Dict[str, int]
                                    ) -> None:
    y_positions = [lineage_map[lineage] for lineage in lineage_map]
    total_values = [mutation_counts.get(lineage, 0) for lineage in lineage_map]
    de_novo_values = [de_novo_counts.get(lineage, 0) for lineage in lineage_map]
    ancestor_values = [t - d for t, d in zip(total_values, de_novo_values)]

    ax.barh(y_positions, de_novo_values, color='red', alpha=0.8, height=0.4,
            edgecolor='black', label='De Novo')
    ax.barh(y_positions, ancestor_values, left=de_novo_values, color='gray',
            alpha=0.6, height=0.4, edgecolor='black', label='Ancestral')

    # Add labels
    for y, t in zip(y_positions, total_values):
        ax.text(t + 0.5, y, str(t), va='center', fontsize=20)

    ax.set_xlabel("Mutation Count", fontsize=28, fontweight='bold')
    ax.set_ylabel("Phage Lineage", fontsize=28, fontweight='bold')
    ax.set_title("De Novo vs. Ancestral Mutations", fontsize=24, fontweight='bold')

    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(max(1, max(total_values) // 5)))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    ax.legend(fontsize=20)


def plot_mutations(df: pd.DataFrame,
                   genes: List[Dict[str, Union[str, int]]],
                   ancestor_phage: str,
                   genome_length: int,
                   output_path: str,
                   keep_lineage_name: bool = False
                   ) -> None:
    """Main function to plot the mutations along the genome with histograms and a heatmap."""

    lineages = list(df['Phage Lineage'].unique())  # Do NOT include the ancestor

    # Adjust y-axis spacing to align all elements properly
    host_bacteria_groups = df.groupby(df['Host Bacteria'].str.strip())["Phage Lineage"].unique()

    # Initialize lineage map with extra spacing
    lineage_map = {}

    y_position = 0
    extra_space = 1  # Extra spacing after each Host Bacteria group

    for host, lineages_in_host in host_bacteria_groups.items():
        for lineage in reversed(lineages_in_host):
            lineage_map[lineage] = y_position
            y_position += 0.5
        y_position += extra_space  # Extra spacing

    first_lineage_per_host = {}
    for host, lineages_in_host in host_bacteria_groups.items():
        if len(lineages_in_host) > 0:
            first_lineage_per_host[lineages_in_host[0]] = host

    gene_y = max(lineage_map.values()) + 1

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
    plot_phage_mutations(ax, df, lineage_map, first_lineage_per_host, genome_length)

    # 🔹 Ensure Histograms are Drawn Correctly
    plot_stacked_mutation_histogram(ax_hist, lineage_map, mutation_counts, de_novo_counts)

    # ✅ Plot the reference genome and gene map last to ensure they are visible
    plot_gene_map(ax, genes, gene_y)

    # # Add legend for gene functions in top-left corner
    # legend_handles = [
    #     mpatches.Patch(color=color, label=label)
    #     for label, color in FUNCTION_COLORS.items()
    # ]
    #
    # ax.legend(
    #     handles=legend_handles,
    #     loc='upper right',
    #     bbox_to_anchor=(1.0, 1.15),  # ⬅️ fine-tune vertical and horizontal position
    #     ncol=1,  # ⬅️ single column (can change to 2+)
    #     fontsize=14,
    #     title="Gene Functions",
    #     title_fontsize=16,
    #     frameon=True
    # )

    # ✅ Re-set y-axis to ensure visibility
    ax.set_ylim(-0.5, gene_y + 1.5)
    ax.set_yticks([])

    print("Gene Y:", gene_y)
    # print("Ancestor Y:", ancestor_y)
    print("Lineage Map:", lineage_map)

    # Extract only used functions from the gene annotations
    used_functions = set(gene["function"].strip().lower() for gene in genes)

    # Filter FUNCTION_COLORS to include only those used
    filtered_function_colors = {func: color for func, color in FUNCTION_COLORS.items() if func in used_functions}

    # Create legend handles for functional colors
    legend_handles = [
        mpatches.Patch(facecolor=color, edgecolor='black', label=func.title())
        for func, color in filtered_function_colors.items()
    ]

    # Place legend outside the figure (right side)
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        bbox_to_anchor=(0.78, 1.02),  # ⬅️ carefully positioned above the right subplot
        fontsize=14,
        title="Gene Functions",
        title_fontsize=16,
        frameon=True,
        ncol=1,  # split legend into two columns to make it more compact
    )

    # Optional: give more space for the legend outside
    fig.subplots_adjust(left=0.1, right=0.95, top=0.90, bottom=0.05, wspace=0.1)

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
        genome_length = get_genome_length(genbank_file)

        # Plot mutations
        plot_mutations(df, genes, f'{ancestor_phage}_reference', genome_length, output_file)

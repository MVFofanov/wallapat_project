from collections import defaultdict
import glob
import os
import re
from typing import Dict, List, Tuple, Union

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


def extract_lineage_info(df: pd.DataFrame, ancestor_phage: str, experiment_version: int = 1) -> pd.DataFrame:
    """Extracts phage lineage and host bacteria name depending on experiment version."""

    if experiment_version == 1:
        # Original approach
        df['Ancestor Phage'] = df['combined_id'].apply(lambda x: x.split('.')[0])
        df['Phage Lineage'] = df['combined_id'].apply(lambda x: '.'.join(x.split('.')[1:])[:-4])

        suffix_regexp = r'P\d+L\d+$'
        df['Host Bacteria'] = df['Phage Lineage'].apply(lambda x: re.sub(suffix_regexp, '', x))

    elif experiment_version == 2:
        # New format: *.consensus.PX_HOST_LY
        df['Phage Lineage'] = df['combined_id'].str.replace(r'\.fasta$', '', regex=True)

        # Host is extracted from the last part (e.g., P8_Y1-16_L3 → Y1-16)
        host_extraction_regex = r'\.consensus\.P\d+_(.*?)_L\d+'
        df['Host Bacteria'] = df['Phage Lineage'].str.extract(host_extraction_regex)[0]

        # Optional: if needed
        df['Ancestor Phage'] = ancestor_phage

    else:
        raise ValueError("Unsupported experiment_version. Use 1 or 2.")

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


def plot_gene_map(ax: matplotlib.axes.Axes,
                  genes: List[Dict[str, Union[str, int]]],
                  gene_y: float,
                  mutation_positions: List[int],
                  label_all: bool = True) -> None:
    """Plots genes along the genome, optionally labeling only those overlapping mutations."""
    print("Plotting gene map at y=", gene_y)  # Debug

    for gene in genes:
        gene_start, gene_end = gene["start"], gene["end"]
        gene_mid = (gene_start + gene_end) / 2

        function = gene.get("function", "unknown function").strip().lower()
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

        # Only label if label_all is True or gene overlaps a mutation
        has_mutation = any(gene_start <= pos <= gene_end for pos in mutation_positions)
        if label_all or has_mutation:
            ax.text(
                gene_mid, gene_y + 0.5,
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
            ax.text(-400, y_pos, first_lineage_per_host[lineage], va='center', fontsize=28, ha='right')

    # Ensure consistent font size for x and y axis ticks
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)


def plot_stacked_mutation_histogram(ax: matplotlib.axes.Axes,
                                    lineage_map: Dict[str, int],
                                    mutation_counts: Dict[str, int],
                                    de_novo_counts: Dict[str, int]
                                    ) -> None:
    """
    Plots a stacked horizontal bar chart of ancestral and de novo mutations.
    Now with ancestral on the bottom and de novo on top.
    """
    y_positions = [lineage_map[lineage] for lineage in lineage_map]
    total_values = [mutation_counts.get(lineage, 0) for lineage in lineage_map]
    de_novo_values = [de_novo_counts.get(lineage, 0) for lineage in lineage_map]
    ancestor_values = [t - d for t, d in zip(total_values, de_novo_values)]

    # First draw ancestral (gray), then de novo (red) on top
    ax.barh(y_positions, ancestor_values, color='gray', alpha=0.6, height=0.4,
            edgecolor='black', label='Ancestral')
    ax.barh(y_positions, de_novo_values, left=ancestor_values, color='red',
            alpha=0.8, height=0.4, edgecolor='black', label='De Novo')

    # Add text labels showing total mutation counts
    for y, t in zip(y_positions, total_values):
        ax.text(t + 0.5, y, str(t), va='center', fontsize=20)

    ax.set_xlabel("Mutation Count", fontsize=28, fontweight='bold')
    ax.set_title("Ancestral vs. De Novo Mutations", fontsize=24, fontweight='bold')

    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(max(1, max(total_values) // 5)))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    ax.legend(fontsize=20)


def prepare_lineage_mapping(df: pd.DataFrame, keep_lineage_name: bool = False) -> Tuple[Dict[str, float], Dict[str, str], Dict[str, str], float]:
    """
    Prepares the lineage_map, first_lineage_per_host, and optionally lineage_labels for plotting.
    Returns:
        lineage_map, first_lineage_per_host, lineage_labels, gene_y_position
    """
    host_bacteria_groups = df.groupby(df['Host Bacteria'].str.strip())["Phage Lineage"].unique()

    lineage_map = {}
    y_position = 0
    extra_space = 1  # Extra spacing after each Host Bacteria group

    for host, lineages_in_host in host_bacteria_groups.items():
        for lineage in reversed(lineages_in_host):
            lineage_map[lineage] = y_position
            y_position += 0.5
        y_position += extra_space

    first_lineage_per_host = {}
    for host, lineages_in_host in host_bacteria_groups.items():
        if len(lineages_in_host) > 0:
            first_lineage_per_host[lineages_in_host[0]] = host

    suffix_regexp = r'P\d+L\d+$'
    lineage_labels = {}

    for lineage in lineage_map:
        if keep_lineage_name:
            lineage_labels[lineage] = lineage
        else:
            lineage_labels[lineage] = re.sub(suffix_regexp, '', lineage)

    gene_y = max(lineage_map.values()) + 1
    return lineage_map, first_lineage_per_host, lineage_labels, gene_y


def compute_mutation_counts(df: pd.DataFrame) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Computes total and de novo mutation counts per lineage.
    Returns:
        mutation_counts, de_novo_counts
    """
    mutation_counts = df.groupby('Phage Lineage')['POS'].count().to_dict()
    de_novo_counts = df[df["mutation_type"] == "de_novo"].groupby("Phage Lineage")["POS"].count().to_dict()
    return mutation_counts, de_novo_counts


def get_used_function_colors(genes: List[Dict[str, Union[str, int]]], function_colors: Dict[str, str]) -> Dict[str, str]:
    """
    Filters the FUNCTION_COLORS to only include those used in the GenBank annotations.
    """
    used_functions = set(gene["function"].strip().lower() for gene in genes)
    return {func: color for func, color in function_colors.items() if func in used_functions}


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

    lineage_map, first_lineage_per_host, lineage_labels, gene_y = prepare_lineage_mapping(df, keep_lineage_name)

    suffix_regexp = r'P\d+L\d+$'
    lineage_labels = {}

    for lineage in lineage_map:
        if keep_lineage_name:
            lineage_labels[lineage] = lineage
        else:
            lineage_labels[lineage] = re.sub(suffix_regexp, '', lineage)

    df['Lineage Order'] = df['Phage Lineage'].map(lineage_map)

    mutation_counts, de_novo_counts = compute_mutation_counts(df)

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

    ax.set_xlabel("Genome position", fontsize=28,
                  fontweight='bold')  # if you have this, or just place below ax.xaxis settings
    ax.set_ylabel("Bacterial host and phage lineages", fontsize=28,
                  fontweight='bold', labelpad=150)
    ax.tick_params(axis='y', which='major', pad=0)  # Reduce distance from ticks to axis

    # 🔹 Ensure Histograms are Drawn Correctly
    plot_stacked_mutation_histogram(ax_hist, lineage_map, mutation_counts, de_novo_counts)

    # ✅ Plot the reference genome and gene map last to ensure they are visible
    # plot_gene_map(ax, genes, gene_y)

    mutation_positions = df["POS"].tolist()
    label_all_genes = len(genes) < 20

    plot_gene_map(ax, genes, gene_y, mutation_positions, label_all=label_all_genes)

    # ✅ Re-set y-axis to ensure visibility
    ax.set_ylim(-0.5, gene_y + 1.5)
    ax.set_yticks([])

    print("Gene Y:", gene_y)
    # print("Ancestor Y:", ancestor_y)
    print("Lineage Map:", lineage_map)

    # Extract only used functions from the gene annotations
    filtered_function_colors = get_used_function_colors(genes, FUNCTION_COLORS)

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
        fontsize=24,
        title="Gene Functions",
        title_fontsize=28,
        frameon=True,
        ncol=1,  # split legend into two columns to make it more compact
    )

    # Optional: give more space for the legend outside
    fig.subplots_adjust(left=0.1, right=0.95, top=0.90, bottom=0.05, wspace=0.1)

    # ✅ Set a title for the entire figure
    ancestor_name = ancestor_phage.rstrip("_reference")
    fig.suptitle(f"{ancestor_name} descendants mutation analysis in different bacterial hosts", fontsize=40, fontweight='bold', y=1.15)

    # 🔹 **Ensure output directory exists and save**
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.savefig(output_path.replace(".png", ".svg"), bbox_inches='tight')  # Save as SVG
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
        df = extract_lineage_info(df, ancestor_phage=ancestor_phage, experiment_version=1)

        # Merge infectivity index
        df = merge_infectivity_data(df, index_df)

        # Load gene annotations
        genes = parse_genbank(genbank_file)
        genome_length = get_genome_length(genbank_file)

        # Plot mutations
        plot_mutations(df, genes, f'{ancestor_phage}_reference', genome_length, output_file)

    # 2025-05-08 new run analysis

    new_fasta_dir = f"{wd}/genomes_new_run_2025_05"
    new_gbk_dir = f"{wd}/genomes_new_run_2025_05_annotation"
    infectivity_dir = f"{wd}/infectivity_index"

    file_path = f"{wd}/20250509_REF_MUT_analysis.xlsx"

    # Get FASTA filenames (e.g., MN988483.P3_L1.consensus.fasta)
    fasta_files = glob.glob(os.path.join(new_fasta_dir, "*.fasta"))
    phages = []
    phage_to_ncbi = {}

    for f in fasta_files:
        base = os.path.basename(f).replace(".consensus.fasta", "")
        ncbi_id, phage_name = base.split(".", 1)
        phages.append(phage_name)
        phage_to_ncbi[phage_name] = ncbi_id

    # Group by pairs like P3_P4 → [P3_L1, P4_L1]

    pair_to_phages = defaultdict(list)
    for p in phages:
        group = "_".join(p.split("_")[:1])  # Use P3, P4 etc
        pair_to_phages[group].append(p)

    index_df_dict = {}  # Cache index tables for each pair

    for group_prefix in ["P3-P4", "P5-P6", "P7-P8"]:
        try:
            index_path = os.path.join(infectivity_dir, f"infectivity_index_{group_prefix}.csv")
            index_df_dict[group_prefix] = load_infectivity_data(index_path)
        except FileNotFoundError:
            print(f"❌ Could not find infectivity file for {group_prefix}")

    for phage in phages:
        ncbi = phage_to_ncbi[phage]
        genbank_file = f"{new_gbk_dir}/{ncbi}.{phage}.gbk"
        sheet_name = f"{phage.replace('_', '')}_filter"
        output_file = f"{figures}/{phage}_phage_mutations_with_genes_de_novo.png"

        group = "-".join(phage[:2] + str(int(phage.split("_")[0][1:]) + 1))  # crude P3_L1 → P3-P4 etc

        index_df = None
        for g, df in index_df_dict.items():
            if phage.startswith(g.split("-")[0]) or phage.startswith(g.split("-")[1]):
                index_df = df
                break
        if index_df is None:
            print(f"⚠️ Skipping {phage} – no infectivity index found")
            continue

        try:
            df = load_mutation_data(file_path, sheet_name)
        except ValueError:
            print(f"⚠️ Sheet {sheet_name} not found in Excel")
            continue

        # df = extract_lineage_info(df, ancestor_phage=phage)
        extract_lineage_info(df, ancestor_phage=phage, experiment_version=2)
        df = merge_infectivity_data(df, index_df)
        genes = parse_genbank(genbank_file)
        genome_length = get_genome_length(genbank_file)

        plot_mutations(df, genes, f'{phage}_reference', genome_length, output_file)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
import urllib.request
import ssl
from matplotlib.lines import Line2D

# Helper function to get the MANE-select transcript
def get_mane_select(transcript_ids):
    mane_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz"
    
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE
    
    with urllib.request.urlopen(mane_url, context=ssl_context) as response:
        mane_df = pd.read_csv(response, sep="\t", compression="gzip")
    
    id_map = {tid.split(".")[0]: tid for tid in transcript_ids}
    mane_df["Ensembl_nuc_clean"] = mane_df["Ensembl_nuc"].str.split(".").str[0]
    
    mane_select = mane_df[mane_df["MANE_status"] == "MANE Select"]["Ensembl_nuc_clean"].values
    
    return {original: (clean in mane_select) for clean, original in id_map.items()}

# Helper function to generate the plot 
def plot_extsp(variant, gene, plot_df, mane_results):
    plot_df = plot_df[(plot_df["Variant"] == variant) & (plot_df["Gene"] == gene)]

    tissue_order = sorted(plot_df["Tissue"].unique())
    transcript_order = sorted(plot_df["Transcript_id"].unique())

    tissue_idx = {t: i for i, t in enumerate(tissue_order)}
    transcript_idx = {tr: i for i, tr in enumerate(transcript_order)}

    x = plot_df["Tissue"].map(tissue_idx)
    y = plot_df["Transcript_id"].map(transcript_idx)
    c = plot_df["exTSP"]

    fig, ax = plt.subplots(figsize=(20, 5))

    norm = mcolors.Normalize(vmin=c.min(), vmax=np.percentile(c, 99))
    sc = ax.scatter(x, y, c=c, cmap=plt.cm.Blues, norm=norm, s=80, edgecolors="steelblue", linewidths=0.5)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("exTSP", fontsize=11)

    ax.set_xticks(range(len(tissue_order)))
    ax.set_xticklabels(tissue_order, rotation=90, fontsize=7)
    ax.set_yticks(range(len(transcript_order)))

    ax.set_yticklabels(transcript_order, fontsize=9)
    for label in ax.get_yticklabels():
        tr = label.get_text()
        if mane_results.get(tr, False):
            label.set_fontweight("bold")

    ax.set_xlabel("GTEx tissues", fontsize=12, fontweight="bold")
    ax.set_ylabel("Transcript", fontsize=12, fontweight="bold")
    ax.set_title(f"{gene} | {variant}\nexTSP across tissues (color = exTSP)", fontsize=13, fontweight="bold")

    ax.set_xlim(-0.5, len(tissue_order) - 0.5)
    ax.set_ylim(-0.5, len(transcript_order) - 0.5)
    ax.grid(True, axis="x", linestyle="--", alpha=0.3)

    # Legend note for MANE select
    legend_elements = [Line2D([0], [0], marker='', color='none', label='Bold = MANE Select transcript')]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.7,
              prop={"weight": "bold", "size": 9})

    plt.tight_layout()
    plt.show()
    return fig, ax

# Load the Cardiomyopathy data
file_path = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/exTSP_file_expression_ensemble_annotation_Cardiomyopathy.tsv"
df = pd.read_csv(file_path, sep='\t')

# Get mane results
mane_results = get_mane_select(df["Transcript_id"].unique())

# Filter for VUS
sub_df = df[df['ClinVar_annotation'] == "Uncertain significance"]

# Then only select those with a high enough IsoPath score to not be LB/B
sub_df = sub_df[sub_df['IsoPath'] > 0.0999]

# Target tissues for Cardiomyopathy
target_tissues = ["Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Coronary",
"Artery_Aorta", "Muscle_skeletal"]

# Filter for all variant/isoform pairs where the max value in relevanttissues is 3x the max in non-relevant tissues
# Store all of these in a dataframe called "hits"
threshold = 3.0
target_mask = sub_df["Tissue"].isin(target_tissues)
target_max = (
    sub_df[target_mask]
    .groupby(["Variant", "Transcript_id"])["exTSP"]
    .max()
    .rename("target_max")
)
rest_max = (
    sub_df[~target_mask]
    .groupby(["Variant", "Transcript_id"])["exTSP"]
    .max()
    .rename("rest_max")
)
means_df = pd.concat([target_max, rest_max], axis=1)
means_df["ratio"] = means_df["target_max"] / means_df["rest_max"]
hitting_combos = means_df[means_df["ratio"] >= threshold].index
hits = sub_df.set_index(["Variant", "Transcript_id"]).loc[hitting_combos].reset_index()

# Filter for variants that have at least five isoforms for plot
sub_hits = hits.groupby("Variant").filter(lambda x: x["Transcript_id"].nunique() >= 5)

# If you want to list out all genes and variants to try different combinations
# print(sub_hits[["Gene", "Variant"]].drop_duplicates())

# Variant for TNNT2
fig, ax = plot_extsp("1_201361314_C/A", "TNNT2", sub_df, mane_results)
# Save as PDF
fig.savefig('outputs/fig3a_TNN2.pdf', bbox_inches='tight')

# Variant for SCN5A
fig, ax = plot_extsp("3_38551307_C/T", "SCN5A", sub_df, mane_results)
# Save as PDF
fig.savefig('outputs/fig3a_SCN5A.pdf', bbox_inches='tight')

# Variant for TPM1
fig, ax = plot_extsp("15_63042855_A/G", "TPM1", sub_df, mane_results)
# Save as PDF
fig.savefig('outputs/fig3a_TPM1.pdf', bbox_inches='tight')

# Variant for MYBPC3
fig, ax = plot_extsp("11_47338588_C/A", "MYBPC3", sub_df, mane_results)
# Save as PDF
fig.savefig('outputs/fig3a_MYBPC3.pdf', bbox_inches='tight')

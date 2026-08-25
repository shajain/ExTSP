# %%
"""Combined Figure 3 — Fig3A (TPM1 scatter) stacked over Fig3B (gene-wise
mini heatmaps), sharing a single GTEx-tissue x-axis.

Fig3A shows transcripts x tissues for one variant (TPM1) as a coloured scatter.
Fig3B shows variant/transcript x tissues per disease as pcolormesh heatmaps.
Both use the same tissue ordering, so they can be laid out one over the other
with a single, shared x-axis (labelled only on the bottom-most panel).
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import urllib.request
import ssl
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP, get_tripletSets_with_exTSP_old

# parse_known_args so the script still runs under interactive/# %% execution
parser = argparse.ArgumentParser(description="Combined Fig3A (TPM1) over Fig3B")
parser.add_argument("--color-code", action="store_true",
                    help="Colour-code disease labels and relevant tissue ticks")
args, _ = parser.parse_known_args()
color_code = args.color_code

use_old_data = True

# Font size for the (shared) GTEx tissue x-axis tick labels — much larger than
# the original per-figure scripts (which used 7 / 9.5).
TISSUE_FONTSIZE = 18

relevant_tissues = {
    "CM":  "Heart_Left_Ventricle",
    "PCD": "Testis",
    "ASD": "Brain_Frontal_Cortex_BA9",
}
colors = {'CM': '#0077B6', 'PCD': '#8261C9', 'ASD': '#D95F02'}
# Map each relevant tissue to its disease colour (for x-axis tick highlighting).
tissue_colors = {relevant_tissues[d]: c for d, c in colors.items()}

# ── Helper functions (shared) ─────────────────────────────────────────────────

def get_mane_select(transcript_ids):
    """Return {transcript_id: is_MANE_Select} for the given transcript ids."""
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


def load_vus(disease, old_data=use_old_data, annotation="Uncertain significance"):
    if old_data:
        df = get_tripletSets_with_exTSP_old(disease, "all")
    else:
        df = get_tripletSets_with_exTSP(disease, "all")
    df["Relevant_Tissue"] = df["Tissue"] == relevant_tissues[disease]
    # For ASD, ClinVar_annotation holds case/control status ("ASD"/"Control")
    # rather than a ClinVar significance, so annotation is passed as "ASD".
    return df[(df['ClinVar_annotation'] == annotation) & (df['IsoPath'] > 0.0999)]


def build_matrix(sub_df, variant_transcript_pairs, all_tissues):
    rows = []
    for variant, transcript in variant_transcript_pairs:
        mask = (sub_df["Variant"] == variant) & (sub_df["Transcript_id"] == transcript)
        vdf  = sub_df[mask][["Tissue", "exTSP"]].groupby("Tissue")["exTSP"].max()
        rows.append(vdf.reindex(all_tissues, fill_value=0.0))
    return np.log1p(pd.DataFrame(rows, columns=all_tissues).values)


# %%
# ── Fig3A data: TPM1 scatter ──────────────────────────────────────────────────

TPM1_VARIANT = "15_63042855_A/G"
TPM1_GENE    = "TPM1"

if use_old_data:
    df_cm = get_tripletSets_with_exTSP_old("CM", "VUS")
else:
    df_cm = get_tripletSets_with_exTSP("CM", "VUS")

mane_results = get_mane_select(df_cm["Transcript_id"].unique())

# Drop likely-benign (LB/B) by IsoPath, matching fig_3a.py
sub_df_3a = df_cm[df_cm['IsoPath'] > 0.0999]
plot_df_3a = sub_df_3a[(sub_df_3a["Variant"] == TPM1_VARIANT) &
                       (sub_df_3a["Gene"] == TPM1_GENE)]

transcript_order = sorted(plot_df_3a["Transcript_id"].unique())
transcript_idx   = {tr: i for i, tr in enumerate(transcript_order)}

# %%
# ── Fig3B data: per-disease heatmaps ──────────────────────────────────────────

# Hard-coded selections (locked in).
hardcoded_variants = {
    "CM": (
        "TPM1",
        [("15_63064114_G/C", "ENST00000403994.9"),
         ("15_63064129_G/C", "ENST00000403994.9"),
         ("15_63042852_T/A", "ENST00000267996.11")],
        ["15_63064114_G/C", "15_63064129_G/C"],
    ),
    "PCD": (
        "DNAAF3",
        [("19_55161656_A/T", "ENST00000524407.7"),
         ("19_55161151_C/G", "ENST00000524407.7"),
         ("19_55165914_G/C", "ENST00000534170.5")],
        ["19_55161656_A/T", "19_55161151_C/G"],
    ),
    "ASD": (
        "SCN2A",
        [("2_165296077_A/G", "ENST00000283256.10"),
         ("2_165313721_G/A", "ENST00000283256.10"),
         ("2_165296068_A/G", "ENST00000631182.3")],
        ["2_165296077_A/G", "2_165313721_G/A"],
    ),
}

diseases = [
    ("CM",  "Uncertain significance"),
    ("PCD", "Uncertain significance"),
    ("ASD", "ASD"),
]

# First pass: load each disease's data (so we can build a shared tissue list).
disease_data = []
for disease, annotation in diseases:
    print(f"\nProcessing {disease}…")
    sub_df = load_vus(disease, annotation=annotation)
    gene, vt_pairs, hit_variants = hardcoded_variants[disease]
    disease_data.append({
        "disease": disease, "sub_df": sub_df, "gene": gene,
        "vt_pairs": vt_pairs, "hit_variants": hit_variants,
    })

# ── Shared tissue ordering ────────────────────────────────────────────────────
# Union of every tissue appearing in Fig3A (TPM1) and Fig3B (all diseases), so
# the scatter and the heatmaps sit on exactly the same x-axis.
tissue_set = set(plot_df_3a["Tissue"].unique())
for d in disease_data:
    tissue_set |= set(d["sub_df"]["Tissue"].unique())
all_tissues  = sorted(tissue_set)
n_tissues    = len(all_tissues)
tissue_idx   = {t: i for i, t in enumerate(all_tissues)}

# Second pass: build heatmap matrices against the shared tissue list.
panels = []
for d in disease_data:
    print(f"\n{d['disease']}: selected gene {d['gene']}")
    for v, t in d["vt_pairs"]:
        tag = 'HIT' if v in d["hit_variants"] else 'NEG-HIT'
        print(f"    [{tag}]  {v}  ->  {t}")
    panels.append({
        "disease"      : d["disease"],
        "gene"         : d["gene"],
        "vt_pairs"     : d["vt_pairs"],
        "hit_variants" : d["hit_variants"],
        "matrix"       : build_matrix(d["sub_df"], d["vt_pairs"], all_tissues),
    })

# Shared colour scale for the heatmaps (log(1 + exTSP)).
all_vals      = np.concatenate([p["matrix"].ravel() for p in panels])
positive_vals = all_vals[all_vals > 0]
vmax_hm = float(np.percentile(positive_vals, 99)) if positive_vals.size else 0.35
norm_hm = mcolors.Normalize(vmin=0, vmax=vmax_hm)

# %%
# ── Combined plot: Fig3A on top, Fig3B (x3) below, shared x-axis ──────────────

n_transcripts = len(transcript_order)

fig_width  = max(24, n_tissues * 0.45)
# Row heights: scatter panel scales with its transcript count; each heatmap is
# 3 variant rows tall. Kept tight so transcripts/cells sit close together.
scatter_h = max(2.4, n_transcripts * 0.28)
heatmap_h = 1.6
fig_height = scatter_h + 3 * heatmap_h + 6.0

fig = plt.figure(figsize=(fig_width, fig_height))

gs = gridspec.GridSpec(
    4, 2,
    height_ratios=[scatter_h, heatmap_h, heatmap_h, heatmap_h],
    width_ratios=[1, 0.015],
    hspace=0.4,
    wspace=0.02,
    left=0.14, right=0.90, top=0.94, bottom=0.26,
)

ax_scatter = fig.add_subplot(gs[0, 0])
ax_heat    = [fig.add_subplot(gs[i, 0], sharex=ax_scatter) for i in (1, 2, 3)]
cbar_ax_a  = fig.add_subplot(gs[0, 1])
cbar_ax_b  = fig.add_subplot(gs[1:, 1])

# ── Fig3A scatter (top) ───────────────────────────────────────────────────────
# Offset x by +0.5 so points sit at the pcolormesh cell centres used below.
x = plot_df_3a["Tissue"].map(tissue_idx) + 0.5
y = plot_df_3a["Transcript_id"].map(transcript_idx)
c = plot_df_3a["exTSP"]

norm_a = mcolors.Normalize(vmin=c.min(), vmax=np.percentile(c, 99))
sc = ax_scatter.scatter(x, y, c=c, cmap=plt.cm.Blues, norm=norm_a,
                        s=90, edgecolors="steelblue", linewidths=0.5)

ax_scatter.set_yticks(range(n_transcripts))
ax_scatter.set_yticklabels(transcript_order, fontsize=15)
for label in ax_scatter.get_yticklabels():
    if mane_results.get(label.get_text(), False):
        label.set_fontweight("bold")
ax_scatter.set_ylabel("Transcript", fontsize=13, fontweight="bold")
ax_scatter.set_ylim(-0.5, n_transcripts - 0.5)
ax_scatter.grid(True, axis="x", linestyle="--", alpha=0.3)

ax_scatter.set_title(f"{TPM1_GENE} | {TPM1_VARIANT}  —  exTSP across tissues",
                     fontsize=15, fontweight="bold")

# Panel label "A"
ax_scatter.text(-0.11, 1.02, "A", transform=ax_scatter.transAxes,
                fontsize=34, fontweight="bold", ha="right", va="bottom")

# Hide scatter x tick labels (shared axis is labelled only at the bottom).
plt.setp(ax_scatter.get_xticklabels(), visible=False)

legend_elements = [Line2D([0], [0], marker='', color='none',
                          label='Bold = MANE Select transcript')]
ax_scatter.legend(handles=legend_elements, loc="upper right", fontsize=10,
                  framealpha=0.7, prop={"weight": "bold", "size": 10})

cbar_a = fig.colorbar(sc, cax=cbar_ax_a)
cbar_a.set_label("exTSP", fontsize=19, rotation=270, labelpad=24)

# ── Fig3B heatmaps (bottom three) ─────────────────────────────────────────────
for idx, (ax, panel) in enumerate(zip(ax_heat, panels)):
    matrix       = panel["matrix"]
    vt_pairs     = panel["vt_pairs"]
    gene         = panel["gene"]
    n_rows       = len(vt_pairs)

    ax.pcolormesh(matrix, cmap="Blues", norm=norm_hm,
                  edgecolors="white", linewidths=0.4)

    # Panel label "B" on the first (top) heatmap.
    if idx == 0:
        ax.text(-0.11, 1.05, "B", transform=ax.transAxes,
                fontsize=34, fontweight="bold", ha="right", va="bottom")

    is_bottom = (idx == len(ax_heat) - 1)
    ax.set_xticks(np.arange(n_tissues) + 0.5)
    if is_bottom:
        ax.set_xticklabels(
            [t.replace("_", " ") for t in all_tissues],
            rotation=90, fontsize=TISSUE_FONTSIZE, ha="center"
        )
        if color_code:
            for tick_label, tissue in zip(ax.get_xticklabels(), all_tissues):
                if tissue in tissue_colors:
                    tick_label.set_color(tissue_colors[tissue])
        ax.set_xlabel("GTEx tissues", fontsize=22, fontweight="bold")
    else:
        plt.setp(ax.get_xticklabels(), visible=False)

    # Y-axis labels (variant bold + transcript italic)
    ax.set_yticks([])
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    for i, (variant, transcript) in enumerate(vt_pairs):
        yc = i + 0.5
        ax.text(-0.01, yc - 0.18, variant,
                transform=trans, ha='right', va='center',
                fontsize=15, fontweight='bold', clip_on=False)
        ax.text(-0.01, yc + 0.18, transcript,
                transform=trans, ha='right', va='center',
                fontsize=13, color='#777777', style='italic', clip_on=False)

    ax.invert_yaxis()

    # Dashed divider between hits and neg-hit (after row 2)
    if n_rows == 3:
        ax.axhline(y=2, color="gray", linewidth=1.0, linestyle="--", alpha=0.6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_edgecolor("black")

    disease_color = colors[panel["disease"]] if color_code else "black"
    ax.text(0.5, 1.04, f'{panel["disease"]}: {gene}',
            transform=ax.transAxes, fontsize=14, fontweight="bold",
            ha="center", va="bottom", color=disease_color,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=disease_color, linewidth=1.5))

# Lock the shared x-range to the pcolormesh extent.
ax_scatter.set_xlim(0, n_tissues)

# Shared heatmap colourbar
cbar_b = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm_hm, cmap="Blues"),
    cax=cbar_ax_b,
)
cbar_b.set_label("log(1 + exTSP)", fontsize=19, rotation=270, labelpad=24)

fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"
output_path = fig_dir + "fig3_combined.pdf"
plt.savefig(output_path, bbox_inches="tight")
print(f"\nSaved to: {output_path}")
plt.show()

# %%
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.transforms import blended_transform_factory
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP, get_tripletSets_with_exTSP_old

# parse_known_args so the script still runs under interactive/# %% execution
parser = argparse.ArgumentParser(description="Fig3B gene-wise mini heatmaps")
parser.add_argument("--color-code", action="store_true",
                    help="Colour-code disease labels and relevant tissue ticks")
args, _ = parser.parse_known_args()
color_code = args.color_code

relevant_tissues = {
    "CM":  "Heart_Left_Ventricle",
    "PCD": "Testis",
    "ASD": "Brain_Frontal_Cortex_BA9",
}
colors = {'CM': '#0077B6', 'PCD': '#8261C9', 'ASD': '#D95F02'}
# Map each relevant tissue to its disease colour (for x-axis tick highlighting).
tissue_colors = {relevant_tissues[d]: c for d, c in colors.items()}
use_old_data = True

# ── Helper functions ────────────────────────────────────────────────────────

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
# ── Per-disease data preparation ─────────────────────────────────────────────

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

panels      = []
all_tissues = None

for disease, annotation in diseases:
    print(f"\nProcessing {disease}…")
    sub_df = load_vus(disease, annotation=annotation)

    if all_tissues is None:
        all_tissues = sorted(sub_df["Tissue"].unique())

    gene, vt_pairs, hit_variants = hardcoded_variants[disease]

    print(f"  Selected gene : {gene}")
    print("  (variant, transcript) pairs:")
    for v, t in vt_pairs:
        print(f"    [{'HIT' if v in hit_variants else 'NEG-HIT'}]  {v}  ->  {t}")

    matrix = build_matrix(sub_df, vt_pairs, all_tissues)
    panels.append({
        "disease"      : disease,
        "gene"         : gene,
        "vt_pairs"     : vt_pairs,
        "hit_variants" : hit_variants,
        "matrix"       : matrix,
    })

# ── Shared colour scale ───────────────────────────────────────────────────────

all_vals      = np.concatenate([p["matrix"].ravel() for p in panels])
positive_vals = all_vals[all_vals > 0]
vmax = float(np.percentile(positive_vals, 99)) if positive_vals.size else 0.35
norm = mcolors.Normalize(vmin=0, vmax=vmax)

# ── Plot ──────────────────────────────────────────────────────────────────────

def _wrap_tissue(t):
    """Split a tissue name into two balanced lines at the nearest character midpoint."""
    words = t.split("_")
    if len(words) <= 1:
        return t
    total = sum(len(w) for w in words)
    half = total / 2
    cumul, best_i, best_diff = 0, 1, float('inf')
    for i, w in enumerate(words[:-1]):
        cumul += len(w)
        diff = abs(cumul - half)
        if diff < best_diff:
            best_diff, best_i = diff, i + 1
    return " ".join(words[:best_i]) + "\n" + " ".join(words[best_i:])

n_tissues    = len(all_tissues)
n_variants   = 3
fig_width    = max(22, n_tissues * 0.32)
panel_height = n_variants * 1.2
fig_height   = panel_height * 3 + 3.5

fig = plt.figure(figsize=(fig_width, fig_height))

gs = gridspec.GridSpec(
    3, 2,
    width_ratios=[1, 0.015],
    hspace=0.35,
    wspace=0.02,
    left=0.12, right=0.92, top=0.90, bottom=0.30,
)

axes    = [fig.add_subplot(gs[i, 0]) for i in range(3)]
cbar_ax = fig.add_subplot(gs[:, 1])

for idx, (ax, panel) in enumerate(zip(axes, panels)):
    matrix        = panel["matrix"]
    vt_pairs      = panel["vt_pairs"]
    hit_variants  = panel["hit_variants"]
    gene          = panel["gene"]

    n_rows = len(vt_pairs)
    ax.pcolormesh(matrix, cmap="Blues", norm=norm,
                  edgecolors="white", linewidths=0.4)

    # X-axis — only the bottom panel gets labels
    ax.set_xticks(np.arange(n_tissues) + 0.5)
    if idx == 2:
        ax.set_xticklabels(
            [_wrap_tissue(t) for t in all_tissues],
            rotation=90, fontsize=9.5, ha="center", fontweight='bold'
        )
        # Colour-code the relevant tissue ticks with their disease colours.
        if color_code:
            for tick_label, tissue in zip(ax.get_xticklabels(), all_tissues):
                if tissue in tissue_colors:
                    tick_label.set_color(tissue_colors[tissue])
        ax.set_xlabel("GTEx tissues", fontsize=13, fontweight="bold")
    else:
        ax.set_xticklabels([])

    # Y-axis labels (variant bold + transcript italic)
    ax.set_yticks([])
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    for i, (variant, transcript) in enumerate(vt_pairs):
        yc = i + 0.5
        ax.text(-0.01, yc - 0.18, variant,
                transform=trans, ha='right', va='center',
                fontsize=11, fontweight='bold', clip_on=False)
        ax.text(-0.01, yc + 0.18, transcript,
                transform=trans, ha='right', va='center',
                fontsize=9.5, color='#777777', style='italic', clip_on=False)

    ax.invert_yaxis()

    # Dashed divider between hits and neg-hit (after row 2)
    if n_rows == 3:
        ax.axhline(y=2, color="gray", linewidth=1.0, linestyle="--", alpha=0.6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_edgecolor("black")

    # Centered label box combining disease + gene, colour-coded per disease
    disease_color = colors[panel["disease"]] if color_code else "black"
    ax.text(0.5, 1.04, f'{panel["disease"]}: {gene}',
            transform=ax.transAxes, fontsize=14, fontweight="bold",
            ha="center", va="bottom", color=disease_color,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=disease_color, linewidth=1.5))

# Shared colourbar
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap="Blues"),
    cax=cbar_ax
)
cbar.set_label("log(1 + exTSP)", fontsize=12, rotation=270, labelpad=14)

fig.suptitle(
    "Fig3B — Gene-wise mini heatmaps (top3 VUS per gene; transcript chosen by mutation's max exTSP)\n"
    "showcasing tissue-specific VUS prioritization",
    fontsize=16, fontweight="bold", y=0.97
)

fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"
output_path = fig_dir + "fig3b.pdf"
plt.savefig(output_path, bbox_inches="tight")
print(f"\nSaved to: {output_path}")
plt.show()

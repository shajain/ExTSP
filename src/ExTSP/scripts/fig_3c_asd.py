# %%

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
import warnings
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP
from ExTSP.isoformSelection.IsoformSelection_Variant import exTSP_selected_isoform

# Shut up the depreciation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

data_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/data/brainspan/"
brainspan_path = data_dir + "ASD_brainspan_by_period_no_P13.tsv"

brainspan_df = pd.read_csv(brainspan_path, sep="\t")
brainspan_df = brainspan_df[brainspan_df['Region'] == 'Cortex']
brainspan_df

# %%
fetal = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
adult = ['P8', 'P9', 'P10', 'P11', 'P12', 'P13']
fetal_df = brainspan_df[(brainspan_df["Developmental_Period"].isin(fetal)) & (brainspan_df["Status"] == "ASD")]
adult_df = brainspan_df[(brainspan_df["Developmental_Period"].isin(adult)) & (brainspan_df["Status"] == "ASD")]


# %%

def valid_labels(case=1, ASD=False):
    '''
    Just get whatever labels fit the case at hand
    '''
    path_labels, benign_labels = [], []
    if ASD:
        assert case < 4
        path_labels = ['ASD']
        benign_labels = ['Control']
    else:
        if case >= 4:
            path_labels = ['Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic']
            benign_labels = ['Benign', 'Likely benign', 'Benign/Likely benign']
        else:
            path_labels = ['Pathogenic']
            benign_labels = ['Benign']
    return path_labels, benign_labels

def select_tissue(df, case=1, ASD=False):
    '''
    Input: The dataframe for the targe disease only

    Let tmax(i) denote the tissue where isoform i has the highest expression.
    In terms of Ptse it can be defined as, tmax(i) = argmaxtPtse(i, t).

    We define S(V, t), the tissue impact score of tissue t, with respect to a
    variant set V , as the number of variant-isoform pairs (v, i) in V × I for
    which t is the tissue with the highest expression of isoform i
    '''

    # Figure out which pathogenic labels we need
    path_labels = valid_labels(case, ASD)[0]
    # Only include pathogenic variants in tissue selection
    df = df[df["ClinVar_annotation"].isin(path_labels)]

    # Select the variants each tissue variant pair with the highest tissue impact probability
    df_subsampled = df.loc[df.groupby(['Transcript_id', 'Variant'])['ptse'].idxmax()]

    # Return the most-impacted tissue
    # print(df_subsampled['Tissue'].value_counts())
    return df_subsampled['Tissue'].value_counts().index[0]

def select_isoform_brainspan(df):
    '''
    Select a single isoform for each variant within each developmental period.

    For every ['Gene', 'Variant', 'Transcript_id', 'Developmental_Period']
    combination, this reduces it to a single
    ['Gene', 'Variant', 'Developmental_Period'] combination by keeping the
    isoform (Transcript_id) whose row has the highest "exTSP" value.
    '''
    return df.loc[df.groupby(["Gene", "Variant", "Developmental_Period"])["exTSP"].idxmax()]

def select_isoform(df, brainspan=False, agg="avg"):
    '''
    Select a single isoform per variant for the exTSP/Tissue metric.

    For brainspan this is a two-phase selection: first select one isoform per
    variant within each developmental period, then aggregate those selections
    over developmental period (avg/max) and select one isoform per variant.
    For GTEx it defers to exTSP_selected_isoform.
    '''
    if not brainspan:
        return exTSP_selected_isoform(df)

    # Phase 1: within each developmental period, select one isoform per variant
    df = pd.concat([
        select_isoform_brainspan(period_df)
        for _, period_df in df.groupby("Developmental_Period")
    ])

    # Phase 2: aggregate the selected isoforms over developmental period,
    # then select one isoform per (Gene, Variant)
    if agg == "avg":
        df = df.groupby(["Gene", "Variant", "Transcript_id"], as_index=False).mean(numeric_only=True)
    else:
        df = df.groupby(["Gene", "Variant", "Transcript_id"], as_index=False).max(numeric_only=True)
    return df.loc[df.groupby(["Gene", "Variant"])["exTSP"].idxmax()]

def proportion_pathogenic_VUS(df, mode="IsoPath", brainspan=False):
    '''
    Compute the proportion of pathogenic VUS for the given metric.

    Expects df to already be ASD-filtered. For mode="Tissue", df must already
    be isoform-selected (see select_isoform); for mode="IsoPath", df should be
    the full (unselected) dataframe.
    '''
    assert mode in ["IsoPath", "Tissue"]

    total_vus = len(df[["Gene", "Variant"]].drop_duplicates())

    if mode == "IsoPath":
        df = df[df['IsoPath'] > 0.0999]
        path_vus = len(df[["Gene", "Variant"]].drop_duplicates())
        return path_vus / total_vus, path_vus, total_vus
    else:
        if brainspan:
            total_vus = len(df)
            df = df[df["exTSP"] > 0.0999]
            path_vus = len(df)
            return path_vus / total_vus, path_vus, total_vus
        else:
            # Select the relevant tissue
            tissue = df[df["Relevant_Tissue"]]["Tissue"].value_counts().index[0]

            df = df[df["Tissue"] == tissue]
            df = df[df["exTSP"] > 0.0999]
            path_vus = len(df[["Gene", "Variant"]].drop_duplicates())
            return path_vus / total_vus, path_vus, total_vus

def plot_fig3c(results):
    """
    Plot Fig3C-style grouped bar chart from a results dictionary.

    Parameters:
        results: dict of the form:
            {
                "Disease Name": {"Tissue": number, "IsoPath": number},
                ...
            }
    """
    diseases = list(results.keys())
    isopath_vals   = [results[d]["IsoPath"][0] * 100 for d in diseases]
    extsp_vals     = [results[d]["Tissue"][0] * 100 for d in diseases]
    isopath_counts = [(results[d]["IsoPath"][1], results[d]["IsoPath"][2]) for d in diseases]
    extsp_counts   = [(results[d]["Tissue"][1], results[d]["Tissue"][2]) for d in diseases]

    x = np.arange(len(diseases))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 7.5))

    bars1 = ax.bar(x - width/2, isopath_vals, width, color="#4472C4", label=r"IsoPath")
    bars2 = ax.bar(x + width/2, extsp_vals,   width, color="#ED7D31", label=r"exTSP")

    for bar, val, (num, den) in zip(bars1, isopath_vals, isopath_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f"{val:.3g}%\n({num}/{den})", ha="center", va="bottom", fontsize=14)

    for bar, val, (num, den) in zip(bars2, extsp_vals, extsp_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f"{val:.3g}%\n({num}/{den})", ha="center", va="bottom", fontsize=14)

    ax.set_ylim(0, 70)
    ax.set_ylabel("VUS re-classified as pathogenic (%)", fontsize=20, y=0.45)
    ax.set_xticks(x)
    ax.set_xticklabels(diseases, rotation=30, ha="right", fontsize=18)
    ax.set_yticks(range(0, 71, 10))
    ax.tick_params(axis="y", labelsize=18)
    ax.legend(fontsize=16, loc="upper right", bbox_to_anchor=(1.0, 1.0))

    ax.set_title(
        "Proportion of VUS reclassified as pathogenic",
        fontsize=20, fontweight="bold", pad=20
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()
    return fig, ax

def main():

    brainspan_df = pd.read_csv(brainspan_path, sep="\t")

    fetal = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
    adult = ['P8', 'P9', 'P10', 'P11', 'P12', 'P13']

    gtex_df = get_tripletSets_with_exTSP("ASD", "case")
    gtex_df['ClinVar_annotation'] = 'ASD'
    gtex_df['Relevant_Tissue'] = gtex_df['Tissue'] == 'Brain_Frontal_Cortex_BA9'

    data_dir = {
        "GTEx": gtex_df,
        "Brainspan Prenatal": brainspan_df[brainspan_df["Developmental_Period"].isin(fetal)],
        "Brainspan Postnatal": brainspan_df[brainspan_df["Developmental_Period"].isin(adult)]
    }

    fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"

    for agg in ["avg", "max"]:
        results = {}
        for disease in data_dir:
            brainspan = "Brainspan" in disease
            status_col = "Status" if brainspan else "ClinVar_annotation"
            df = data_dir[disease]
            df = df[df[status_col] == "ASD"]

            # All isoform selection happens here, before proportion_pathogenic_VUS
            selected_df = select_isoform(df, brainspan=brainspan, agg=agg)

            disease_results = {
                "IsoPath": proportion_pathogenic_VUS(df, mode="IsoPath", brainspan=brainspan),
                "Tissue": proportion_pathogenic_VUS(selected_df, mode="Tissue", brainspan=brainspan),
            }
            results[disease] = disease_results

        fig, ax = plot_fig3c(results)
        fig.savefig(fig_dir + f'fig3c_asd_{agg}.pdf', bbox_inches='tight')


if __name__ == "__main__":
    main()

# %%

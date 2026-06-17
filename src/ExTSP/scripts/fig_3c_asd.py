# %%

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
import warnings
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP

# Shut up the depreciation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

data_dir = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/"
brainspan_path = data_dir + "exTSP_file_expression_ensemble_annotation_ASD_Satt_Fu_BrainSpan.tsv"

brainspan_df = pd.read_csv(brainspan_path, sep="\t")
fetal = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
adult = ['P8', 'P9', 'P10', 'P11', 'P12', 'P13']
fetal_df = brainspan_df[(brainspan_df["Tissue"].isin(fetal)) & (brainspan_df["ClinVar_annotation"] == "ASD")]
adult_df = brainspan_df[(brainspan_df["Tissue"].isin(adult)) & (brainspan_df["ClinVar_annotation"] == "ASD")]


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

def select_isoform(df):
    '''
    For each variant we consider the isoform with the highest value of Ptse(i, t) as the most
    important isoform for the variant. We refer to this isoform as Ptse-select isoform.

    Returns a reduced version of the dataframe with the selected isoform for each variant
    '''
    # Step 1: Calculate the max ptse score for each Variant-Transcript_id combination
    transcript_scores = df.groupby(['Variant', 'Transcript_id'])['ptse'].max().reset_index()
    transcript_scores.columns = ['Variant', 'Transcript_id', 'max_ptse']

    # Step 2: For each Variant, find the Transcript_id with the highest max_ptse
    # Shuffle first so idxmax breaks ties randomly (pandas 3.0 compat: avoid apply with group key)
    shuffled = transcript_scores.sample(frac=1)
    idx = shuffled.groupby('Variant')['max_ptse'].idxmax()
    best_transcripts = shuffled.loc[idx, ['Variant', 'Transcript_id']].reset_index(drop=True)

    # Step 3: Filter the original dataframe to keep only rows with the selected Variant-Transcript_id pairs
    df_reduced = df.merge(best_transcripts, on=['Variant', 'Transcript_id'], how='inner')
    return df_reduced

def proportion_pathogenic_VUS(df, mode="IsoPath", brainspan=False, agg="avg"):
    assert mode in ["IsoPath", "Tissue"]
    assert agg in ["avg", "max"]

    df = df[df['ClinVar_annotation'] == "ASD"]
    total_vus = len(df[["Gene", "Variant"]].drop_duplicates())

    if mode == "IsoPath":
        df = df[df['IsoPath'] > 0.0999]
        path_vus = len(df[["Gene", "Variant"]].drop_duplicates())
        return path_vus / total_vus, path_vus, total_vus
    else:
        # df = select_isoform(df)
        if brainspan:
            ## Aggregate over tissues (avg or max)
            if agg == "avg":
                df = df.groupby(["Gene", "Variant", "Transcript_id"], as_index=False).mean(numeric_only=True)
                df = df.loc[df.groupby(["Gene", "Variant"])["exTSP"].idxmax()]
            else:
                df = df.groupby(["Gene", "Variant", "Transcript_id"], as_index=False).max(numeric_only=True)
                df = df.loc[df.groupby(["Gene", "Variant"])["exTSP"].idxmax()]
            total_vus = len(df) #BEcause we've already removed redundant tissues
            df = df[df["exTSP"] > 0.0999]
            path_vus = len(df)
            return path_vus / total_vus, path_vus, total_vus
        else:
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

    fig, ax = plt.subplots(figsize=(12, 7))

    bars1 = ax.bar(x - width/2, isopath_vals, width, color="#4472C4", label=r"IsoPath $\geq$ 0.0999 (relevant rows)")
    bars2 = ax.bar(x + width/2, extsp_vals,   width, color="#ED7D31", label=r"exTSP $\geq$ 0.0999 (relevant rows)")

    for bar, val, (num, den) in zip(bars1, isopath_vals, isopath_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f"{val:.3g}%\n({num}/{den})", ha="center", va="bottom", fontsize=8.5)

    for bar, val, (num, den) in zip(bars2, extsp_vals, extsp_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f"{val:.3g}%\n({num}/{den})", ha="center", va="bottom", fontsize=8.5)

    ax.set_ylim(0, 70)
    ax.set_ylabel("Pathogenic VUS (%)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(diseases, rotation=30, ha="right", fontsize=11)
    ax.set_yticks(range(0, 71, 10))
    ax.legend(fontsize=10, loc="upper left")

    ax.set_title(
        "Fig3C —Proportion of VUS reclassified as pathogenic across ASD Datasets \n"
        "IsoPath vs exTSP (Relevant_Tissue and max over transcripts)",
        fontsize=13, fontweight="bold"
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
        "Brainspan Prenatal": brainspan_df[brainspan_df["Tissue"].isin(fetal)],
        "Brainspan Postnatal": brainspan_df[brainspan_df["Tissue"].isin(adult)]
    }

    fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"

    for agg in ["avg", "max"]:
        results = {}
        for disease in data_dir:
            df = data_dir[disease]
            brainspan = "Brainspan" in disease
            disease_results = {}
            for mode in ["IsoPath", "Tissue"]:
                ratio, num, den = proportion_pathogenic_VUS(df, mode=mode, brainspan=brainspan, agg=agg)
                disease_results[mode] = (ratio, num, den)
            results[disease] = disease_results

        fig, ax = plot_fig3c(results)
        fig.savefig(fig_dir + f'fig3c_asd_{agg}.pdf', bbox_inches='tight')


if __name__ == "__main__":
    main()

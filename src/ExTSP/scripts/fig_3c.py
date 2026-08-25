import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP
from ExTSP.isoformSelection.IsoformSelection_Variant import exTSP_selected_isoform
import warnings

# Shut up the depreciation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

data_dir = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/"
asd_path = data_dir + "exTSP_file_expression_ensemble_annotation_ASD_Satt_Fu_GTEx.tsv"
cm_path = data_dir + "exTSP_file_expression_ensemble_annotation_Cardiomyopathy.tsv"
kd_path = data_dir + "exTSP_file_expression_ensemble_annotation_KD.tsv"
pcd_path = data_dir + "exTSP_file_expression_ensemble_annotation_PCD.tsv"
ph_path = data_dir + "exTSP_file_expression_ensemble_annotation_PH.tsv"

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

def proportion_pathogenic_VUS(df, mode="IsoPath"):
    assert mode in ["IsoPath", "Tissue select", "Relevant Tissue"]

    print(f"Running Isoform Selection")
    df = exTSP_selected_isoform(df)
    if mode == "Tissue select":
        tissue = select_tissue(df)
    df = df[df['ClinVar_annotation'] == "Uncertain significance"]
    total_vus = len(df[["Gene", "Variant"]].drop_duplicates())
    
    if mode == "IsoPath":
        df = df[df['IsoPath'] > 0.0999]
        path_vus = len(df[["Gene", "Variant"]].drop_duplicates())
        return path_vus / total_vus, path_vus, total_vus
    else:
        if mode == "Relevant Tissue":
            # filter for relevant tissue
            tissue = df[df["Relevant_Tissue"]]["Tissue"].value_counts().index[0]
        df = df[df["Tissue"] == tissue]
        total_vus = len(df) #BEcause we've already removed redundant tissues
        df = df[df["exTSP"] > 0.0999]
        path_vus = len(df)
        return path_vus / total_vus, path_vus, total_vus

def plot_fig3c(results):
    """
    Plot Fig3C-style grouped bar chart from a results dictionary.

    Parameters:
        results: dict of the form:
            {
                "Disease Name": {"Relevant Tissue": number, "IsoPath": number},
                ...
            }
    """
    diseases = list(results.keys())
    isopath_vals   = [results[d]["IsoPath"][0] * 100 for d in diseases]
    extsp_vals     = [results[d]["Relevant Tissue"][0] * 100 for d in diseases]
    isopath_counts = [(results[d]["IsoPath"][1], results[d]["IsoPath"][2]) for d in diseases]
    extsp_counts   = [(results[d]["Relevant Tissue"][1], results[d]["Relevant Tissue"][2]) for d in diseases]

    x = np.arange(len(diseases))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 7))

    bars1 = ax.bar(x - width/2, isopath_vals, width, color="#4472C4", label=r"IsoPath")
    bars2 = ax.bar(x + width/2, extsp_vals,   width, color="#ED7D31", label=r"exTSP")

    for bar, val, (num, den) in zip(bars1, isopath_vals, isopath_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                f"{val:.3g}%\n({num}/{den})", ha="center", va="bottom", fontsize=13)

    for bar, val, (num, den) in zip(bars2, extsp_vals, extsp_counts):
        ax.text(bar.get_x() + bar.get_width()/2 + 0.05, bar.get_height() + 1.5,
                f"{val:.3g}%\n({num}/{den})", ha="center", va="bottom", fontsize=13)

    ax.set_ylim(0, 60)
    ax.set_ylabel("Proportion of VUS reclassified as pathogenic (%)", fontsize=20, y=0.50)
    ax.set_xticks(x)
    ax.set_xticklabels(diseases, rotation=30, ha="right", fontsize=23)
    ax.set_yticks(range(0, 51, 10))
    ax.tick_params(axis="y", labelsize=20)
    ax.tick_params(axis="x", length=0)
    ax.legend(fontsize=16, loc="upper right")

    ax.set_title(
        "Proportion of VUS reclassified as pathogenic across diseases",
        fontsize=20, fontweight="bold", pad=20
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig("fig3c.png", dpi=150)
    plt.show()
    return fig, ax

def main():
    # Old (unnecesary)
    # file_dir = {
    #     "CM": cm_path,
    #     "PCD": pcd_path,
    #     "PH": ph_path,
    #     "KD": kd_path
    # }

    # New implementation
    diseases = ['CM', 'PCD', 'PH', 'PKD']
    relevant_tissues = {
        'CM':  'Heart_Left_Ventricle',
        'PCD': 'Testis',
        'PH':  'Lung',
        'PKD': 'Kidney_Medulla',
    }

    results = {}
    for disease in diseases:
        # df = pd.read_csv(file_dir[disease], sep="\t") #Old
        print(f"Loading {disease} data")
        df = get_tripletSets_with_exTSP(disease, 'all')
        df['Relevant_Tissue'] = df['Tissue'] == relevant_tissues[disease]
        disease_results = {}
        for mode in ["IsoPath", "Tissue select", "Relevant Tissue"]:
            disease_results[mode] = proportion_pathogenic_VUS(df, mode=mode)
        results[disease] = disease_results

    fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"

    fig, ax = plot_fig3c(results)
    fig.savefig(fig_dir + 'fig3c.pdf', bbox_inches='tight')


if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import pickle
import sys
import argparse
from ExTSP.commonFunctions import summaryStats, filterPathOrCases, clinSig_bool
from ExTSP.config import DATA_DIR, RESULTS_DIR, FIGURES_DIR, NUM_CPU, NBOOT
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP

PERFORM_STATISTICAL_TESTS = True



def ComputeTissueEnrichment(df):
    Tissues = df['Tissue'].unique()
    #df.Transcript_id.nunique()
    print(df.shape)
    # Given each Variant-Transcript pair, get the tissue with the maximum ExTSP value. 
    max_ix = df.groupby(['Variant', 'Transcript_id'])['ptse'].idxmax()
    df_max_Tissue = df.loc[max_ix]
    DF_tissue  = []
    i = 0
    VTPair_counts = []
    for _, group_df in df_max_Tissue.groupby('Gene'):
        #print("\n\nGene:", group_df['Gene'].iloc[0])
        # For each tissue count the number of variant-transcripts where it has max ExTSP value
        Tissue_counts = group_df['Tissue'].value_counts().to_frame()
        #print(Tissue_counts)
        # Set the 0 counts for tissues that don't get picked
        remainingTissues = set(Tissues) - set(Tissue_counts.index.to_list())
        for tissue in remainingTissues:
            Tissue_counts.loc[tissue] = 0
        total_count = Tissue_counts.sum().to_numpy()[0].item()
        Tissue_counts["proportion"] = Tissue_counts / total_count
        VTPair_counts.append(total_count)
        DF_tissue.append(Tissue_counts)
        if i == 0:
            df_tissue = Tissue_counts
        else:   
            df_tissue[["count", "proportion"]] += Tissue_counts
        i += 1
    df_tissue["proportion"]= df_tissue['proportion'] / len(DF_tissue)
    df_tissue["proportion_og"]= df_tissue['count'] / df_tissue['count'].sum()
    return df_tissue

def ComputeTissueEnrichment(df):
    Tissues = df['Tissue'].unique()
    #df.Transcript_id.nunique()
    print(df.shape)
    # Given each Variant-Transcript pair, get the tissue with the maximum ExTSP value. 
    max_ix = df.groupby(['Variant', 'Transcript_id'])['ptse'].idxmax()
    df_max_Tissue = df.loc[max_ix]
    DF_tissue  = []
    i = 0
    VTPair_counts = []
    for _, group_df in df_max_Tissue.groupby('Gene'):
        #print("\n\nGene:", group_df['Gene'].iloc[0])
        # For each tissue count the number of variant-transcripts where it has max ExTSP value
        Tissue_counts = group_df['Tissue'].value_counts().to_frame()
        #print(Tissue_counts)
        # Set the 0 counts for tissues that don't get picked
        remainingTissues = set(Tissues) - set(Tissue_counts.index.to_list())
        for tissue in remainingTissues:
            Tissue_counts.loc[tissue] = 0
        total_count = Tissue_counts.sum().to_numpy()[0].item()
        Tissue_counts["proportion"] = Tissue_counts / total_count
        VTPair_counts.append(total_count)
        DF_tissue.append(Tissue_counts)
        if i == 0:
            df_tissue = Tissue_counts
        else:   
            df_tissue[["count", "proportion"]] += Tissue_counts
        i += 1
    df_tissue["proportion"]= df_tissue['proportion'] / len(DF_tissue)
    df_tissue["proportion_og"]= df_tissue['count'] / df_tissue['count'].sum()
    return df_tissue


# Null distribution
def genNullDistribution(df_ntarget, key, maxTissue, nboot=20, num_cpu=8):
    print(maxTissue)
    print("nboot:", nboot)
    print("num_cpu:", num_cpu)
    df_VT = df_ntarget[['Variant', 'Transcript_id']].drop_duplicates()
    maxTissueProp = []
    def function(b):
        print("Bootstrap:", b)
        #df_sampled_VT = df_VT.sample(n=df_VT.shape[0], replace=True)
        #df_ntarget_sampled = df_ntarget.merge(df_sampled_VT, on=['Variant', 'Transcript_id'], how='inner')
        idx = np.random.choice(len(df_ntarget),len(df_ntarget),replace=True)
        df_ntarget_sampled = df_ntarget.iloc[idx]
        enrich_ntarget = ComputeTissueEnrichment(df_ntarget_sampled)[key].to_frame().rename(columns={key: f'boot_{b}'})
        return enrich_ntarget
    DFs = Parallel(n_jobs=num_cpu)(delayed(function)(b) for b in range(nboot))
    enrich_ntarget = pd.concat(DFs, axis=1)
    maxTissueProp = [enrich_ntarget.loc[maxTissue, f'boot_{b}'].item() for b in range(nboot)]
    return enrich_ntarget, maxTissueProp



def filterPathOrCases(df, disease):
    if disease == "ASD" and "Status" in df.columns:
        df = df[df.Status == "ASD"]
    else:
        P_idx, PLP_idx = clinSig_bool(clinSig=df.ClinicalSignificance)[0:2]
        df = df[P_idx|PLP_idx]
    return df


def getTargetAndNonTargetDFs_old(disease, Files):
    VIT_df = pd.read_csv(Files[disease]['VIT'], sep='\t')
    #VIT_df = VIT_df[idxPLP(VIT_df.ClinVar_annotation)]
    targetVars_df = pd.read_csv(Files[disease]['targetVars'], sep='\t')
    targetVars_df = filterPathOrCases(targetVars_df, disease)
    nonTargetVars_df = pd.read_csv(Files[disease]['nontargetVars'], sep='\t')
    nonTargetVars_df = filterPathOrCases(nonTargetVars_df, disease)
    target_df = VIT_df[VIT_df['Variant'].isin(targetVars_df['Variant'])]
    nonTarget_df = VIT_df[VIT_df['Variant'].isin(nonTargetVars_df['Variant'])]
    # print("non-target ", disease, nonTarget_df.shape)
    # nonTarget_df = nonTarget_df[~nonTarget_df['Variant'].isin(targetVars_df['Variant'])]
    # print("non-target ", disease, nonTarget_df.shape)
    # nonTarget_df = nonTarget_df[~nonTarget_df['Gene'].isin(targetVars_df['Gene'])]
    # print("non-target ", disease, nonTarget_df.shape)
    return target_df, nonTarget_df

def getTargetAndNonTargetDFs(disease):
    type = "PLP" if disease!="ASD" else "case"
    targetVars_df = get_tripletSets_with_exTSP(disease, type=type)
    nonTargetVars_df = get_tripletSets_with_exTSP('nonTarget', type="PLP")
    return targetVars_df, nonTargetVars_df

def lolipop_plot(df_target, df_ntarget, statSig, key, col, title):
    maxTissue = df_target[key].idxmax()
    fig, ax = plt.subplots(figsize=(10, 2) )
    df_target = df_target.sort_index()
    print(df_target[key].sum())
    df_ntarget = df_ntarget.sort_index()
    print(df_ntarget[key].sum())
    assert np.isclose(df_target[key].sum(), 1.0), "Target proportions do not sum to 1!"
    assert np.isclose(df_ntarget[key].sum(), 1.0), "Non-Target proportions do not sum to 1!"
    markerline, stemlines, baseline = ax.stem(df_target[key] ,label='Target', markerfmt='o', basefmt=" ")
    # Style markers
    markerline.set_markerfacecolor(col)
    markerline.set_markeredgecolor('none')
    markerline.set_markersize(5)
    #markerline.set_marker('o')

    # Style stem lines
    stemlines.set_color('black')
    stemlines.set_linewidth(1)
    #stemlines.set_linestyle('-')
    markerline, stemlines, baseline = ax.stem(df_ntarget[key] ,label='Non-Target', markerfmt='o',linefmt='none', basefmt=" ")
    markerline.set_markerfacecolor("0.7")
    markerline.set_markeredgecolor('none')
    markerline.set_markersize(5)

    if statSig:
        ax.plot(df_target.index.to_list().index(maxTissue), df_target.loc[maxTissue, key].item()+0.07, marker='*', color='black', markersize=5)
    ax.set_ylabel('Proportion of variant-transcript pairs')
    tissues = df_target.index.to_list()
    assert tissues == df_ntarget.index.to_list(), "Tissue lists do not match!"
    ax.set_xticks(range(len(tissues)))
    ax.set_xticklabels(tissues, rotation=90, ha ='center')
    ax.set_ylim(0, 1.05)
    plt.title(title)
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/{title}_tissueEnrichment.pdf", bbox_inches='tight')
    #plt.show()

def extractSDVRows(df):
    df_SV = df[df["SV"]==1]
    return df_SV 


if __name__ == "__main__":
    # if len(sys.argv) > 1:
    #     NBOOT = int(sys.argv[1])
    #     print("Number of bootstraps set to:", NBOOT)
    # if len(sys.argv) > 2:
    #     NUM_CPU = int(sys.argv[2])
    parser = argparse.ArgumentParser(description="")

    # Add keyword arguments
    parser.add_argument("--NBOOT", type=int, help="Number of bootstraps for statistical significance", required=False)
    parser.add_argument("--NUM_CPU", type=int, help="Number of CPUs to use (default: 1).", required=False)

    # Parse the arguments
    args = parser.parse_args()
    if args.NBOOT:
        NBOOT = args.NBOOT
        PERFORM_STATISTICAL_TESTS = True
    if args.NUM_CPU:
        NUM_CPU = args.NUM_CPU
    print("Number of bootstraps set to:", NBOOT)
    print("Number of CPUs set to:", NUM_CPU)
        
    #diseases = ["CM", "PKD", "PH", "PCD", "ASD"]
    diseases = ["PKD", "PH", "PCD", "ASD"]
    #diseases = ["ASD"]
    files = {d: {"VIT": f'{DATA_DIR}/TissueEnrichment/{d}_VITs.tsv', 
                 "targetVars": f'{DATA_DIR}/TissueEnrichment/{d}_targetVars.tsv', 
                 "nontargetVars": f'{DATA_DIR}/TissueEnrichment/{d}_nontargetVars.tsv'} for d in diseases}
    colors= {"CM": '#8B3649', "PKD": '#d95f02', "PH": '#7570b3', "PCD": '#4468A7', "ASD": '#839741'}
    #key = 'proportion'
    key = 'proportion_og'
    for disease in diseases:
        # target_df, nonTarget_df = getTargetAndNonTargetDFs(disease, files)
        target_df, nonTarget_df = getTargetAndNonTargetDFs(disease)
        # if disease == "CM":
        #     print(target_df.shape)
        #     VCEP = ["MYH7","MYBPC3","TNNT2","TPM1","MYL2","MYL3","TNNI3","ACTC1"]
        #     target_df = target_df[target_df["Gene"].isin(VCEP)]
        #     print(target_df.shape)
        #     #nonTarget_df = extractSDVRows(nonTarget_df)
        summaryStats(target_df, disease, target=True, type="SDV")
        summaryStats(nonTarget_df, disease, target=False, type="SDV")
        TE_target = ComputeTissueEnrichment(target_df)
        TE_ntarget = ComputeTissueEnrichment(nonTarget_df)
        pickle.dump( (TE_target, TE_ntarget), 
                    open( f"{RESULTS_DIR}/{disease}_tissueEnrichment.pkl", "wb" ) )
        maxTissue = TE_target[key].idxmax()
        sig = False
        if PERFORM_STATISTICAL_TESTS:
            TE_ntarget_boot, maxTissueProp_in_ntarget = genNullDistribution(nonTarget_df, key, maxTissue, nboot=NBOOT, num_cpu=NUM_CPU)
            sig = TE_target.loc[maxTissue, key].item() >= np.percentile(maxTissueProp_in_ntarget, 95)
            pickle.dump( (sig, maxTissueProp_in_ntarget,TE_ntarget_boot), 
                        open( f"{RESULTS_DIR}/{disease}_nullDist.pkl", "wb" ) )  

        lolipop_plot(TE_target, TE_ntarget, sig, key, colors[disease], title=f"{disease}")

        
        

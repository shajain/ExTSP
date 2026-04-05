import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import pdb
import seaborn as sns
import argparse
from ExTSP.isoformSelection.IsoformSelection_Gene import SUPPORTING_EVIDENCE_THRESHOLD
from ExTSP.commonFunctions import filterPathOrCases, summaryStats
import seaborn as sns

SUPPORTING_EVIDENCE_THRESHOLD = 0.0999


def exTSP_selected_isoform_Tissue(df, Tissue):
    df = df[df['Tissue'] == Tissue]
    return exTSP_selected_isoform(df)

def exTSP_selected_isoform(df, selectionColumn="exTSP"):
    df_exTSP_selected = pd.DataFrame()
    for _, group_df in df.groupby('Variant'):
        #get the row with maximum exTSP
        exTSP_selected_transcript = group_df.loc[group_df[selectionColumn].idxmax()]["Transcript_id"]
        group_df_exTSP_selected = group_df[group_df["Transcript_id"] == exTSP_selected_transcript]
        df_exTSP_selected = pd.concat([df_exTSP_selected, group_df_exTSP_selected])
    return df_exTSP_selected.reset_index(drop=True)

def isoformSelectionVariant(df, bestTissue, relevantTissues, useBestTissueTranscript=True, topk=None, topr=2):
    df_Path = filterPathOrCases(df)
    Tissues = df['Tissue'].unique()
    nonRevTissue_ranks = pd.DataFrame(columns=['Tissue', 'rank'])
    nonRelTissues = sorted(list(set(Tissues) - set(relevantTissues+[bestTissue])))
    nonRevTissue_ranks['Tissue'] = nonRelTissues
    nonRevTissue_ranks['rank'] = 0
    nonRevTissue_ranks['exTSPAvg'] = 0.0
    nonRelcols = ['nonRevTissue'+str(i) + s for i in range(1, topr+1) for s in ['', 'Trannscript', 'ExTSP']] 
    #columns=["Gene", 'Variant', 'bestTissue', 'bestTissueTranscript', 'bestExTsp', 'nonRevTissue1', 'nonRevTissue1Transcript', 'nonRevTissue1ExTsp', 
    #                                          'nonRevTissue2', 'nonRevTissue2Transcript', 'nonRevTissue2ExTsp']
    columns=["Gene", 'Variant', 'bestTissue', 'bestTissueTranscript', 'bestExTSP'] +  nonRelcols
    df_varHeatmap = pd.DataFrame(columns=columns)
    df_varHeatmap1 = df_varHeatmap.copy()
    df_varHeatmap1 = df_varHeatmap.copy()
    DF_remaining = pd.DataFrame()
    DF_best = pd.DataFrame()
    displayVar_df = None
    cnts = 0
    for _, group_df in df_Path.groupby('Gene'):
        group_df_cpy = group_df.copy()
        gene = group_df["Gene"].to_list()[0]
        idx = group_df[(group_df["Tissue"]==bestTissue) & (group_df["exTSP"]>=SUPPORTING_EVIDENCE_THRESHOLD)].groupby('Variant')["exTSP"].idxmax()
        var_df = group_df.loc[idx].sort_values(by='exTSP', ascending=False)
        if topk is not None:
            var_df = var_df[0:topk]
        Variant_list = set(var_df["Variant"].to_list())
        #var_df = var_df[["Variant", "exTSP"]].drop_duplicates().sort_values(by='IsoPath', ascending=False)
        #vars_filtered = var_df["Variant"].unique().tolist()
        group_df = group_df[group_df["Variant"].isin(Variant_list)]
        
       
        if useBestTissueTranscript:
            # For best tissue, select the transcript with highest exTSP. Also pick only those transcripts for other tissues.
            max_ix_bt = group_df[group_df['Tissue'] == bestTissue].groupby(['Variant'])["exTSP"].idxmax()
            group_df_bt = group_df.loc[max_ix_bt]
            group_df_other = group_df[group_df['Tissue'] != bestTissue]
            group_df_other = group_df_other.merge(group_df_bt[['Variant', 'Transcript_id']], on=['Variant', 'Transcript_id'], how='right')
            group_df = pd.concat([group_df_bt, group_df_other])
        else:
            max_ix = group_df.groupby(['Variant', 'Tissue'])["exTSP"].idxmax()
            group_df = group_df.loc[max_ix]
        df_best = group_df[group_df['Tissue'] == bestTissue]
        df_relevant = group_df[group_df['Tissue'].isin(relevantTissues)]
        df_remaining = group_df[~group_df['Tissue'].isin(set(relevantTissues).union({bestTissue}))]
        DF_best = pd.concat([DF_best, df_best])
        DF_remaining = pd.concat([DF_remaining, df_remaining])
        for _, df_var_group in df_remaining.groupby('Variant'):
            df_var_group.sort_values(by='exTSP', ascending=False, inplace=True)
            df_var_group.loc[:, 'rank'] = range(1, len(df_var_group) + 1)
            ranks = df_var_group.sort_values(by='Tissue', inplace=False)['rank'].reset_index(drop=True)
            exTSP = df_var_group.sort_values(by='Tissue', inplace=False)['exTSP'].reset_index(drop=True)
            nonRevTissue_ranks.loc[:, 'rank'] += ranks==1
            nonRevTissue_ranks.loc[:, 'exTSPAvg'] += exTSP
            print(ranks.loc[29],ranks.loc[30])
            cnts+=1
        #nonRevTissue_ranks.loc[:, 'rank'] /= topk
        # var_df = df_best[["Variant", "exTSP"]].sort_values(by='exTSP', ascending=False)[0:topk]
        # Variant_list = var_df["Variant"].to_list().unique()
        # # group_df = group_df.merge(var_df[["Variant"]], on=['Variant'], how='inner')
        df_top_nonRel = df_remaining.groupby(['Variant'], group_keys=False).apply(lambda x: x.sort_values(by="exTSP", ascending=False).head(topr)).reset_index(drop=True)
        #df_top_nonRel = df_top_nonRel[df_top_nonRel['Variant'].isin(Variant_list)]
        for var in Variant_list:
            entriesFromBest = df_best[df_best['Variant'] == var].iloc[0][['Tissue', 'Transcript_id', 'exTSP']].to_list()
            entriesFromNonRel = df_top_nonRel[df_top_nonRel['Variant'] == var].loc[:, ['Tissue', 'Transcript_id', 'exTSP']].to_numpy().flatten().tolist()
            row = [gene, var] + entriesFromBest + entriesFromNonRel
            df_varHeatmap.loc[len(df_varHeatmap)] = row
    nonRevTissue_ranks.loc[:, 'rank'] /= cnts
    nonRevTissue_ranks.loc[:, 'exTSPAvg'] /= cnts
    #global_top_nonRel = nonRevTissue_ranks.sort_values(by='rank').head(topr)["Tissue"].to_list()
    global_top_nonRel = nonRevTissue_ranks.sort_values(by='exTSPAvg', ascending=False).head(topr)["Tissue"].to_list()
    df_top_nonRel1 = DF_remaining[DF_remaining['Tissue'].isin(global_top_nonRel)]
    for idx,row in df_varHeatmap.iterrows():
        gene = row["Gene"]
        var = row["Variant"]
        entriesFromBest = row[['bestTissue', 'bestTissueTranscript', 'bestExTSP']].tolist()
        entriesFromNonRel = df_top_nonRel1[(df_top_nonRel1["Gene"]==gene) & (df_top_nonRel1["Variant"]==var)][['Tissue', 'Transcript_id', 'exTSP']]
        entriesFromNonRel = entriesFromNonRel.set_index('Tissue').reindex(global_top_nonRel).reset_index()
        entriesFromNonRel = entriesFromNonRel.to_numpy().flatten().tolist()
        row = [gene, var] + entriesFromBest + entriesFromNonRel
        df_varHeatmap1.loc[len(df_varHeatmap1)] = row

    df_PATH_MANEselect = df_Path[df_Path['MANE_Select']][["Gene", "Variant", "Transcript_id"]].drop_duplicates()
    df_PATH_MANEselect  = df_PATH_MANEselect.rename(columns={'Transcript_id': 'MANE_Select'})
    df_Path_MANEClinical = df_Path[df_Path['MANE_Plus_Clinical']][["Gene", "Variant", "Transcript_id"]].drop_duplicates()
    df_Path_MANEClinical = df_Path_MANEClinical.rename(columns={'Transcript_id': 'MANE_Plus_Clinical'})
    df_varHeatmap =  df_varHeatmap.merge(df_PATH_MANEselect, on=['Gene', 'Variant'], how='left')
    df_varHeatmap =  df_varHeatmap.merge(df_Path_MANEClinical, on=['Gene', 'Variant'], how='left')
    df_varHeatmap1 =  df_varHeatmap1.merge(df_PATH_MANEselect, on=['Gene', 'Variant'], how='left')
    df_varHeatmap1 =  df_varHeatmap1.merge(df_Path_MANEClinical, on=['Gene', 'Variant'], how='left')
    
    # df_varHeatmap2 = DF_best.pivot_table(index=['Gene', 'Variant'], columns='Tissue', values='exTSP').reset_index()
    # df_varHeatmap2 = df_varHeatmap2.merge(DF_remaining.pivot_table(index=['Gene', 'Variant'], columns='Tissue', values='exTSP').reset_index(), on=['Gene', 'Variant'], how='inner')
    return df_varHeatmap, df_varHeatmap1, displayVar_df


def plot_heatmap(df_HM, topr=2, transcriptIdentical=True):
# Extract exTSP scores for heatmap
    bestTissue = df_HM['bestTissue'].unique()[0]
    nonRelcols = ['nonRevTissue'+str(i) + s for i in range(1, topr+1) for s in ['', 'Trannscript', 'ExTSP']] 
    columns = ['bestExTSP'] + [nonRelcols[2+i*3] for i in range(topr)]
    heatmap_data = df_HM[columns].values

    # Create y-axis labels with Gene, Variant, Transcript info
    y_labels = [f"{row['Gene']}|{row['Variant']}|{row['bestTissueTranscript']}" 
                for _, row in df_HM.iterrows()]


    # Create annotations with tissue names, transcripts, and scores
    annotations = []
    for _, row in df_HM.iterrows():
        if not transcriptIdentical:
             row_annot = [ f"{row['bestTissue']}\n{row['bestTissueTranscript']}\n{row['bestExTSP']:.3f}"] \
             + [ f"{row[nonRelcols[0+i*3]]}\n{row[nonRelcols[1+i*3]]}\n{row[nonRelcols[2+i*3]]:.3f}" for i in range(topr)]
        else:
            val =  f"{row['bestExTSP']:.3f}" if row["bestExTSP"] > SUPPORTING_EVIDENCE_THRESHOLD else ""
            #val =  f"{row['bestExTSP']:.3f}"
            row_annot = [val] + [ f"{row[nonRelcols[2+i*3]]:.3f}" if row[nonRelcols[2+i*3]] > SUPPORTING_EVIDENCE_THRESHOLD else "" for i in range(topr)]
        annotations.append(row_annot)
    if transcriptIdentical:
        Tissues = [bestTissue] + [df_HM[nonRelcols[i]].unique()[0] for i in range(0, len(nonRelcols), 3)]
    else:
        Tissues = [bestTissue] + ['Tissue '+ str(i) for i in range(1,topr)]
    # Create heatmap
    #plt.figure(figsize=(14, len(df_varHeatmap1) * 0.6))
    if transcriptIdentical:
        fig, ax = plt.subplots(figsize=(10, len(df_HM) * 0.3))
    else:
        fig, ax = plt.subplots(figsize=(30, len(df_HM) * 0.6))
    sns.heatmap(heatmap_data, 
                annot=annotations,    # Show tissue, transcript, and score
                fmt='',               # Empty format since we're using custom annotations
                cmap='Blues',         # Color scheme
                cbar_kws={'label': 'exTSP Score'},
                xticklabels= Tissues,
                yticklabels=y_labels)

    plt.xlabel('Tissue Category')
    plt.ylabel('Gene | Variant | Transcript')
    plt.title('exTSP Scores Across Tissues')
    for i, (idx, row) in enumerate(df_HM.iterrows()):
        transcript = row['bestTissueTranscript']
        MANE_SELECT = row['MANE_Select']
        MANE_CLINICAL = row['MANE_Plus_Clinical']
        # Check if transcript is MANE_Select or MANE_Plus_Clinical
        mane_not_present = MANE_SELECT is np.nan and MANE_CLINICAL is np.nan
        not_mane = transcript != MANE_SELECT and transcript != MANE_CLINICAL
        yticklabel = ax.get_yticklabels()[i]
        if mane_not_present:
            yticklabel.set_color('#4468A7' )
        elif not_mane:
            yticklabel.set_color('#950606' )  # MANE transcripts in green
    plt.tight_layout()
    plt.savefig(f"Figures/IsoformSelection/Variant_Heatmap_{disease}1.pdf", bbox_inches='tight')
    #plt.show()



if __name__ == "__main__":

   
    parser = argparse.ArgumentParser(description="")

    # # Add keyword arguments
    parser.add_argument("--disease", type=int, help="Pick disease from CM, KD, PH, PCD, ASD", required=False)
    BestTissue = {"CM": "Heart_Left_Ventricle","KD":"Kidney_Medula","PH":"Liver","PCD":"Lung","ASD":"Brain_Cortex"}    

    # Parse the arguments
    args = parser.parse_args()
    if args.disease:
        disease = args.disease
    else:
        disease = "CM"

    bestTissue = BestTissue[disease]

    print("Disease set to:", disease)
    print("Best tissue set to:", bestTissue)
    
    vitFile = f'data/IsoformSelection/{disease}_VITs.tsv'
    key = 'proportion_og'
    vit_df = pd.read_csv(vitFile, sep='\t')
    summaryStats(vit_df, disease, target=True, type="MDV")
    topk = 2
    topr = 5
    df_varHeatmap, df_varHeatmap1, displayVar_df = isoformSelectionVariant(vit_df, bestTissue=bestTissue, 
                                                                     relevantTissues=[bestTissue], 
                                                                     useBestTissueTranscript=True, 
                                                                     topk=topk, topr=topr)
    #plot_heatmap(df_varHeatmap, topr=topr, transcriptIdentical=False)
    plot_heatmap(df_varHeatmap1, topr=topr)
    
        
    

        #lolipop_plot(TE_target, TE_ntarget, sig, key, colors[disease], title=f"{disease}")

        
        

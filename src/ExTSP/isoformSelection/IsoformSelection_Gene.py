import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import pdb
import seaborn as sns
import argparse
from ExTSP.commonFunctions import filterPathOrCases, summaryStats, getSingleDiseaseVariants
import seaborn as sns
from matplotlib.colors import to_rgba


SUPPORTING_EVIDENCE_THRESHOLD = 0.0999


# def exTSP_selected_isoform(df, Tissue=None):
#     if Tissue is not None:
#         df = df[df['Tissue'] == Tissue]
#     for _, group_df in df.groupby('Gene'):
#         df_IT = group_df.groupby(["Transcript_id", "Tissue"])["exTSP"].mean()
#         df_isoform_sorted = df_IT.groupby("Transcript_id")["exTSP"].max().sort_values(by='exTSP', ascending=False).reset_index()[["Transcript_id"]]
#         Variants = group_df["Variant"].unique().tolist()
#         rank = 1
#         while len(Variants) > 0:

#     return df_isoform


def isoformSelectionGene(df, bestTissue, relevantTissues, useBestTissueTranscript=True, topk=None, topr=2):
    df_Path = filterPathOrCases(df, onlyP=False)
    Tissues = df['Tissue'].unique()
    nonRevTissue_ranks = pd.DataFrame(columns=['Tissue', 'rank'])
    nonRelTissues = sorted(list(set(Tissues) - set(relevantTissues+[bestTissue])))
    nonRevTissue_ranks['Tissue'] = nonRelTissues
    nonRevTissue_ranks['rank'] = 0
    nonRevTissue_ranks["exTSPAvg"] = 0.0
    nonRevTissue_ranks["nSupportingEvi"] = 0
    nonRelcols = ['nonRevTissue'+str(i) + s for i in range(1, topr+1) for s in ['', 'Trannscript', 'ExTSP']] 
    #columns=["Gene", 'Variant', 'bestTissue', 'bestTissueTranscript', 'bestExTsp', 'nonRevTissue1', 'nonRevTissue1Transcript', 'nonRevTissue1ExTsp', 
    #                                          'nonRevTissue2', 'nonRevTissue2Transcript', 'nonRevTissue2ExTsp']
    columns=["Gene", 'bestTissue', 'bestTissueTranscript', 'bestExTSP'] +  nonRelcols
    df_geneHeatmap = pd.DataFrame(columns=columns)
    df_geneHeatmap1 = df_geneHeatmap.copy()
    df_geneHeatmap1 = df_geneHeatmap.copy()
    DF_remaining = pd.DataFrame()
    DF_best = pd.DataFrame()
    cnts = 0
    for _, group_df in df_Path.groupby('Gene'):
        gene = group_df["Gene"].to_list()[0]
        #idx = group_df[(group_df["Tissue"]==bestTissue) & (group_df["exTSP"]>=0.0999)].groupby('Variant')["exTSP"].idxmax()
        #var_df = group_df.loc[idx]
        #var_df = var_df["Tissue", ].to_frame().rename(columns={'Tissue': 'bestTissue'})
        # if topk is not None:
        #     var_df = var_df[0:topk]
        # Variant_list = set(var_df["Variant"].to_list())
        #group_df = group_df[group_df["Variant"].isin(Variant_list)]
        if not group_df[group_df['Tissue'] == bestTissue]["exTSP"].max() > SUPPORTING_EVIDENCE_THRESHOLD:
            continue
        TT_df = group_df[["Tissue", "Transcript_id", "exTSP"]].groupby(['Tissue','Transcript_id'])["exTSP"].mean().reset_index()

        if useBestTissueTranscript:
            # For best tissue, select the transcript with highest exTSP. Also pick only those transcripts for other tissues.
            max_ix_bt = TT_df[TT_df['Tissue'] == bestTissue]["exTSP"].idxmax()
            bestTranscript = TT_df.loc[max_ix_bt]['Transcript_id']
            TT_df = TT_df[TT_df['Transcript_id'] == bestTranscript]
            # group_df_bt = group_df.loc[max_ix_bt]
            # group_df_other = group_df[group_df['Tissue'] != bestTissue]
            # group_df_other = group_df_other.merge(group_df_bt[['Variant', 'Transcript_id']], on=['Variant', 'Transcript_id'], how='right')
            # group_df = pd.concat([group_df_bt, group_df_other])
        else:
            TT_df = TT_df[TT_df.groupby('Tissue')["exTSP"].idxmax()]
            # max_ix = group_df.groupby(['Variant', 'Tissue'])["exTSP"].idxmax()
            # group_df = group_df.loc[max_ix]
        df_best = TT_df[TT_df['Tissue'] == bestTissue]
        df_relevant = TT_df[TT_df['Tissue'].isin(relevantTissues)]
        df_remaining = TT_df[~TT_df['Tissue'].isin(set(relevantTissues).union({bestTissue}))]
        DF_best = pd.concat([DF_best, df_best])
        df2_remaining = df_remaining.copy()
        df2_remaining['Gene'] = gene
        DF_remaining = pd.concat([DF_remaining, df2_remaining])
        df_remaining.sort_values(by='exTSP', ascending=False, inplace=True)
        df_remaining.loc[:, 'rank'] = range(1, len(df_remaining) + 1)
        ranks = df_remaining.sort_values(by='Tissue', inplace=False)['rank'].reset_index(drop=True)
        exTSP = df_remaining.sort_values(by='Tissue', inplace=False)['exTSP'].reset_index(drop=True)
        nonRevTissue_ranks.loc[:, 'rank'] += ranks==1 
        nonRevTissue_ranks.loc[:, 'exTSPAvg'] += exTSP
        nonRevTissue_ranks.loc[:, 'nSupportingEvi'] += exTSP>SUPPORTING_EVIDENCE_THRESHOLD
        cnts+=1
        #nonRevTissue_ranks.loc[:, 'rank'] /= topk
        # var_df = df_best[["Variant", "exTSP"]].sort_values(by='exTSP', ascending=False)[0:topk]
        # Variant_list = var_df["Variant"].to_list().unique()
        # # group_df = group_df.merge(var_df[["Variant"]], on=['Variant'], how='inner')
        
        #df_top_nonRel = df_top_nonRel[df_top_nonRel['Variant'].isin(Variant_list)]
        entriesFromBest = df_best.iloc[0].tolist()
        df_top_nonRel = df_remaining[["Tissue", "Transcript_id", "exTSP"]].sort_values(by="exTSP", ascending=False).head(topr).reset_index(drop=True)
        entriesFromNonRel = df_top_nonRel.to_numpy().flatten().tolist()
        row = [gene] + entriesFromBest + entriesFromNonRel
        df_geneHeatmap.loc[len(df_geneHeatmap)] = row
    nonRevTissue_ranks.loc[:, 'rank'] /= cnts
    nonRevTissue_ranks.loc[:, 'exTSPAvg'] /= cnts
    #global_top_nonRel = nonRevTissue_ranks.sort_values(by='rank', ascending=False).head(topr)["Tissue"].to_list()
    global_top_nonRel = nonRevTissue_ranks.sort_values(by='exTSPAvg', ascending=False).head(topr)["Tissue"].to_list()
    df_top_nonRel1 = DF_remaining[DF_remaining['Tissue'].isin(global_top_nonRel)]
    for idx,row in df_geneHeatmap.iterrows():
        gene = row["Gene"]
        entriesFromBest = row[['bestTissue', 'bestTissueTranscript', 'bestExTSP']].tolist()
        entriesFromNonRel = df_top_nonRel1[df_top_nonRel1["Gene"] == gene][["Tissue", "Transcript_id", "exTSP"]]
        entriesFromNonRel = entriesFromNonRel.set_index('Tissue').reindex(global_top_nonRel).reset_index()
        entriesFromNonRel = entriesFromNonRel.to_numpy().flatten().tolist()
        row = [gene] + entriesFromBest + entriesFromNonRel
        df_geneHeatmap1.loc[len(df_geneHeatmap1)] = row

    df_PATH_MANEselect = df_Path[df_Path['MANE_Select']][["Gene", "Transcript_id"]].drop_duplicates()
    df_PATH_MANEselect  = df_PATH_MANEselect.rename(columns={'Transcript_id': 'MANE_Select'})
    df_Path_MANEClinical = df_Path[df_Path['MANE_Plus_Clinical']][["Gene", "Transcript_id"]].drop_duplicates()
    df_Path_MANEClinical = df_Path_MANEClinical.rename(columns={'Transcript_id': 'MANE_Plus_Clinical'})
    df_geneHeatmap =  df_geneHeatmap.merge(df_PATH_MANEselect, on=['Gene'], how='left')
    df_geneHeatmap =  df_geneHeatmap.merge(df_Path_MANEClinical, on=['Gene'], how='left')
    df_geneHeatmap1 =  df_geneHeatmap1.merge(df_PATH_MANEselect, on=['Gene'], how='left')
    df_geneHeatmap1 =  df_geneHeatmap1.merge(df_Path_MANEClinical, on=['Gene'], how='left')
    
    # df_geneHeatmap2 = DF_best.pivot_table(index=['Gene', 'Variant'], columns='Tissue', values='exTSP').reset_index()
    # df_geneHeatmap2 = df_geneHeatmap2.merge(DF_remaining.pivot_table(index=['Gene', 'Variant'], columns='Tissue', values='exTSP').reset_index(), on=['Gene', 'Variant'], how='inner')
    return df_geneHeatmap, df_geneHeatmap1



def isoformSelectionExample(df, bestTissue):
    #filter to pathogenic variants and best tissue
    df_Path = filterPathOrCases(df, onlyP=False)
    df_Path_bt = df_Path[df_Path['Tissue'] == bestTissue]

    output = {}
    for gene, df_gene in df_Path_bt.groupby('Gene'):
        if gene == "TPM1" or gene == "FHL1" or gene == "CACNA1C" or gene == "ACTA1":
            print(gene)
        #Skip if not enough variants or transcripts
        nVariants = df_gene['Variant'].nunique()
        nTranscripts = df_gene['Transcript_id'].nunique()
        if nTranscripts >= 10 and df_gene["exTSP"].max() > SUPPORTING_EVIDENCE_THRESHOLD:
            print(df_gene[df_gene["exTSP"] == df_gene["exTSP"].max()]["Variant"].to_list())
        if nVariants < 3 or nTranscripts < 2:
            continue
        # Make a list of MANE transcripts. Skip if none present.
        mane_select = df_gene[df_gene["MANE_Select"]]["Transcript_id"].unique().tolist()
        mane_plus_clinical = df_gene[df_gene["MANE_Plus_Clinical"]]["Transcript_id"].unique().tolist()
        mane_transcripts = mane_select + mane_plus_clinical
        if len(mane_transcripts) == 0:
            continue

        #Compute average exTSP per transcript over variants and count the variants supporting each transcript
        varCount_df = df_gene[["Transcript_id", "Variant"]].groupby('Transcript_id')["Variant"].nunique().reset_index()
        varCount_df.rename(columns={'Variant': 'Variant_counts'}, inplace=True)
        varAvg_df = df_gene[["Transcript_id", "exTSP"]].groupby('Transcript_id')["exTSP"].mean().reset_index()
        varAvg_df = varAvg_df.merge(varCount_df, on='Transcript_id', how='left')
        varAvg_df.sort_values(by='exTSP', ascending=False, inplace=True)

        #Skip if average exTSP of best tissue is not greater than supporting evidence threshold
        # if varAvg_df["exTSP"].max() <= SUPPORTING_EVIDENCE_THRESHOLD:
        #     continue

        print("Gene:", gene)
        print("Number of variants:", nVariants)
        print("Number of transcripts:", nTranscripts)   
        print("max average exTSP:", varAvg_df["exTSP"].max())
        print("MANE max average exTSP:", varAvg_df[varAvg_df["Transcript_id"].isin(mane_transcripts)]["exTSP"].max())

        #Skip gene if avg exTSP of non-MANE is not higher than MANE
        max_AvgExTSP_mane = varAvg_df[varAvg_df["Transcript_id"].isin(mane_transcripts)]["exTSP"].max()
        if varAvg_df["exTSP"].max() - max_AvgExTSP_mane <= 0.0:
            continue

        # Find the best display variant with maximum delta between best non-MANE and best MANE exTSP
        deltaMANE = 0.0
        bestVar_df = None
        bestVar = None
        best_max_exTSP = 0.0
        print(df_gene["exTSP"].max())
        for var, df_var in df_gene.groupby('Variant'):
            # variant not present on MANE transcripts skip
            max_exTSP_mane = df_var[(df_var["Transcript_id"].isin(mane_transcripts))]["exTSP"].max()
            if max_exTSP_mane is np.nan:
                continue
            # Skip variant if best exTSP less than supporting evidence threshold
            max_exTSP = df_var["exTSP"].max()
            # if max_exTSP <= SUPPORTING_EVIDENCE_THRESHOLD:
            #     continue
            # Update best display variant if deltaMANE improved
            #if max_exTSP - max_exTSP_mane > deltaMANE:
            if (max_exTSP > max_exTSP_mane) and (max_exTSP > best_max_exTSP):
                deltaMANE = max_exTSP - max_exTSP_mane
                best_max_exTSP = max_exTSP
                bestVar = var
                bestVar_df = df_var[["Transcript_id", "exTSP","MANE_Select","MANE_Plus_Clinical"]].sort_values(by='exTSP', ascending=False)  
        #If no suitable variant found, skip gene
        # if bestVar is None:
        #     continue
        # Annotate varAvg_df with MANE info
        varAvg_df["MANE_Select"] = varAvg_df["Transcript_id"].isin(mane_select)
        varAvg_df["MANE_Plus_Clinical"] = varAvg_df["Transcript_id"].isin(mane_plus_clinical)
        #Create output entry for the gene
        output[gene] = {"varAvg_df": varAvg_df, "bestVar": bestVar, "bestVar_df": bestVar_df}
    return output



def plot_heatmap(df_HM, topr=2, transcriptIdentical=True):
# Extract exTSP scores for heatmap
    bestTissue = df_HM['bestTissue'].unique()[0]
    nonRelcols = ['nonRevTissue'+str(i) + s for i in range(1, topr+1) for s in ['', 'Trannscript', 'ExTSP']] 
    columns = ['bestExTSP'] + [nonRelcols[2+i*3] for i in range(topr)]
    heatmap_data = df_HM[columns].values

    # Create y-axis labels with Gene, Variant, Transcript info
    y_labels = [f"{row['Gene']}|{row['bestTissueTranscript']}" 
                for _, row in df_HM.iterrows()]


    # Create annotations with tissue names, transcripts, and scores
    annotations = []
    for _, row in df_HM.iterrows():
        if not transcriptIdentical:
             row_annot = [ f"{row['bestTissue']}\n{row['bestTissueTranscript']}\n{row['bestExTSP']:.3f}"] \
             + [ f"{row[nonRelcols[0+i*3]]}\n{row[nonRelcols[1+i*3]]}\n{row[nonRelcols[2+i*3]]:.3f}" for i in range(topr)]
        else:
            val =  f"{row['bestExTSP']:.3f}" if row["bestExTSP"] > 0.0999 else ""
            #val =  f"{row['bestExTSP']:.3f}"
            row_annot = [val] + [ f"{row[nonRelcols[2+i*3]]:.3f}" if row[nonRelcols[2+i*3]] > 0.0999 else "" for i in range(topr)]
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
    plt.ylabel('Gene | Transcript')
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
    plt.savefig(f"Figures/IsoformSelection/{disease}/Gene_Heatmap.pdf", bbox_inches='tight')
    #plt.show()

def plot_exampleGeneAndVariants(gene, geneData, disease):
    varAvg_df = geneData["varAvg_df"]
    bestVar = geneData["bestVar"]
    bestVar_df = geneData["bestVar_df"]
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7, 8))
    def loliplop(ax, df, title, yLabel, var=None):
        col =  '#8B3649'
        df = df.sort_values(by='exTSP', ascending=False)
        exTSP_scores = df['exTSP'].values
        transcripts = df['Transcript_id'].values
        mane_select = df["MANE_Select"].values
        mane_clinical = df["MANE_Plus_Clinical"].values
        color_list = ['#4468A7' if mane else col for mane in mane_select | mane_clinical]
        
        # Plot stem for each point with its color
        for i, (score, color) in enumerate(zip(exTSP_scores, color_list)):
            markerline, stemlines, baseline = ax.stem([i], [score], markerfmt='o', basefmt=" ")
            markerline.set_markerfacecolor(color)
            markerline.set_markeredgecolor('none')
            markerline.set_markersize(10)
            stemlines.set_color('black')
            stemlines.set_linewidth(1)
        if var is not None:
            ax.set_title(f'{disease} | {gene} | {var}')
        else:
            ax.set_title(f'{disease} | {gene}')
        ax.set_ylabel(yLabel)
        ax.set_xlabel("Transcripts")
        ax.set_xticks(range(len(transcripts)))
        if var is not None:
            xticklabels = transcripts
        else:
            counts = df['Variant_counts'].values
            xticklabels = [f"{t} ({c})" for (t,c) in zip(transcripts, counts)]
        ax.set_xticklabels(xticklabels, rotation=45, ha='right')
        #ax.set_ylim(0, 1.05)
        #set tight layout for ax  
    loliplop(axs[0], varAvg_df, title=f"{disease} | {gene}", yLabel="Average exTSP")
    if bestVar is not None:
        loliplop(axs[1], bestVar_df, title=f"{disease} | {gene} | {bestVar}", yLabel="exTSP", var=bestVar)  
    plt.tight_layout()
    #How to save the sublots togenther
    plt.savefig(f"Figures/IsoformSelection/{disease}/Example_{gene}_{disease}.pdf", bbox_inches='tight')
    #plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    # # Add keyword arguments
    parser.add_argument("--diseases", type=int, help="Pick disease from CM, KD, PH, PCD, ASD", required=False)
    BestTissue = {"CM": "Heart_Left_Ventricle","KD":"Kidney_Medulla","PH":"Lung","PCD":"Testis","ASD":"Brain_Cerebellar_Hemisphere"}    

    # Parse the arguments
    args = parser.parse_args()
    if args.diseases:
        diseases = args.diseases
    else:
        diseases = ["CM","KD","PH","PCD","ASD"]
    
    for disease in diseases:
        bestTissue = BestTissue[disease]

        print("Disease set to:", disease)
        print("Best tissue set to:", bestTissue)
    
        vitFile = f'data/IsoformSelection/{disease}_VITs.tsv'
        key = 'proportion_og'
        vit_df = pd.read_csv(vitFile, sep='\t')
        #vit_df, _ = getSingleDiseaseVariants(disease)
        summaryStats(vit_df, disease, target=True, type="MDV")
        topk = 3
        topr = 5
        df_varHeatmap, df_varHeatmap1 = isoformSelectionGene(vit_df, bestTissue=bestTissue, 
                                                                     relevantTissues=[bestTissue], 
                                                                     useBestTissueTranscript=True, 
                                                                     topk=topk, topr=topr)
        #plot_heatmap(df_varHeatmap, topr=topr, transcriptIdentical=False)
        plot_heatmap(df_varHeatmap1, topr=topr)
        output = isoformSelectionExample(vit_df, bestTissue)
        for gene in output:
            plot_exampleGeneAndVariants(gene, output[gene], disease)
        
    

        #lolipop_plot(TE_target, TE_ntarget, sig, key, colors[disease], title=f"{disease}")

        
        

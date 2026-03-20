import pandas as pd
import numpy as np

def summaryStats(df, disease, target, type):
    df_GV = df[["Gene", "Variant"]].drop_duplicates()
    num_genes = df_GV['Gene'].nunique()
    num_variants = df_GV['Variant'].nunique()
    if target:
        s = "target"
    else:
        s = "non-target"
    print("Summary statistics for", s, disease)    
    #print("Number of Gene-variant pairs for", disease, ":", num_GV)
    print("Nuber Genes:", num_genes)
    print("Nuber Variants:", num_variants)  
    print("\n")
    return 

def filterPathOrCases(df, onlyP=False):
    if 'ASD' in df.ClinVar_annotation.unique(): 
        df_Path = df.query("ClinVar_annotation in ['ASD']")
    elif onlyP:
        df_Path = df.query("ClinVar_annotation in ['Pathogenic']")
    else:
        df_Path = df.query("ClinVar_annotation in ['Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic']")
    return df_Path  

def clinSig_bool(clinSig, onlyP=False):
    PLP_idx = np.array(['pathogenic' in cs.lower() for cs in clinSig])
    BLB_idx = np.array(['benign' in cs.lower() for cs in clinSig])
    VUS_idx = np.array(['uncertain significance' in cs.lower() for cs in clinSig])  
    likely_idx = np.array(['likely' in cs.lower() for cs in clinSig])
    PLP_idx = PLP_idx & ~BLB_idx & ~VUS_idx
    BLB_idx = BLB_idx & ~PLP_idx & ~VUS_idx
    P_idx = PLP_idx & ~likely_idx
    B_idx = BLB_idx & ~likely_idx
    return P_idx, PLP_idx, VUS_idx, BLB_idx, B_idx

def isCase_ASD(status):
    iscase = np.array([st=="ASD" for st in status])
    return iscase



def getSingleDiseaseVariants(disease):
    diseases = ["CM", "KD", "PH", "PCD", "ASD"]
    Files = {d: {"VIT": f'data/TissueEnrichment/{d}_VITs.tsv', 
                 "targetVars": f'data/TissueEnrichment/{d}_targetVars.tsv', 
                 "nontargetVars": f'data/TissueEnrichment/{d}_nontargetVars.tsv'} for d in diseases}
    VIT_df = pd.read_csv(Files[disease]['VIT'], sep='\t')
    targetVars_df = pd.read_csv(Files[disease]['targetVars'], sep='\t')
    nonTargetVars_df = pd.read_csv(Files[disease]['nontargetVars'], sep='\t')
    target_df = VIT_df[VIT_df['Variant'].isin(targetVars_df['Variant'])]
    nonTarget_df = VIT_df[VIT_df['Variant'].isin(nonTargetVars_df['Variant'])]
    # print("non-target ", disease, nonTarget_df.shape)
    # nonTarget_df = nonTarget_df[~nonTarget_df['Variant'].isin(targetVars_df['Variant'])]
    # print("non-target ", disease, nonTarget_df.shape)
    # nonTarget_df = nonTarget_df[~nonTarget_df['Gene'].isin(targetVars_df['Gene'])]
    # print("non-target ", disease, nonTarget_df.shape)
    return target_df, nonTarget_df
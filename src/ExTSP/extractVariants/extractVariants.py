#from turtle import st
from urllib import response
import pandas as pd
import numpy as np
#import json
#import re
from ExTSP.commonFunctions import clinSig_bool
#import requests
#from config import DISEASE_SEMANTIC_TYPES
from ExTSP.phenotypes.common_functions import find_in_ClinVar_master, find_in_disease_Phenotypes
from ExTSP.extractVariants.usefulFuncs import extractPhenoIDs
from ExTSP.config import Ontologies, DATA_DIR

#UMLS_API_KEY = "a76354ec-df97-4be1-82da-3e3eadd09f46"


Raw_folder = f'{DATA_DIR}/raw'
ClinVar_file = f'{Raw_folder}/missense_variants_ClinVar_Nov_2024.tsv'
# ASDVariantFiles = [f'{Raw_folder}/Antaki_Sebat_SPARK_missence_SupplementaryTable2.xlsx',
#                    f'{Raw_folder}/Fu_Talkowski_Suppl_Table20_ASD_genes_72_2022_preprocessed_Hg38_start_end.xlsx',
#                    f'{Raw_folder}/Satterstrom_TableS1_missense_PTV_variants.xlsx']

filtered_folder = f'{DATA_DIR}/Filtered'
ASDVarFile = f'{filtered_folder}/ASD_Vars.tsv'
extractedVariantsFile = f'{filtered_folder}/extractedVariants.xlsx'



# def onlyKeepDiseaseMedgenIDs(PIDs):
#     medgenInCV = set([id.split(':')[1] for idlist in PIDs for id in idlist if id.split(':')[0].lower() == "medgen"])
#     filtered_medgen_ids = filterMedgenIDs(medgenInCV)
#     exclude_ids = medgenInCV - set(filtered_medgen_ids)
#     print(f"Excluding {len(exclude_ids)} MedGen IDs that do not correspond to diseases: {exclude_ids}")
#     PIDs1 = [[pid for pid in pidList if not (pid.split(':')[0].lower() == "medgen" and pid.split(':')[1] in exclude_ids)] for pidList in PIDs]
#     return PIDs1

def consideredAndMatched(infoCV, infoDisease, ont):
    considered = True if infoCV else False
    matched = True if (infoCV and infoDisease) else False
    assert (not matched) or considered
    if ont.lower() == "medgen":
        isFinding = False
        isAncestors = False
        semantic_type = infoCV.get("semantictype", "")
        #semantic_type_disease = infoDisease.get("semantictype","") 
        hierarchy = infoDisease.get("hierarchy","")
        if type(semantic_type)==str and semantic_type.lower()== "finding":
            isFinding = True
        if type(hierarchy)==str and hierarchy.lower()== "ancestor":
            isAncestors = True
        if isFinding or isAncestors:
            considered = False
            matched = False
    return considered, matched



def phenoInfo2bool(IDs, disease):
    match = sum([len(ids)>0 for ont, ids in IDs.items()])>=1
    if not match:
        assert not match
    total_match_count = 0    
    for ont in IDs.keys():
        if ont.lower() == "omim":
            continue
        ids = IDs[ont]
        inClinVarPheno = [find_in_ClinVar_master(id, ontology=ont) for id in ids]
        inDiseasePheno = [find_in_disease_Phenotypes(id, ontology=ont, disease=disease) for id in ids]
        num_ids_considered = 0
        num_ids_matched = 0
        for infoCV, infoDiease in zip(inClinVarPheno, inDiseasePheno):
            considered, matched = consideredAndMatched(infoCV, infoDiease, ont)
            num_ids_considered += considered
            num_ids_matched +=matched
            total_match_count += matched
        if num_ids_considered >=4:
            # match &= np.floor(num_ids_matched/num_ids_considered) > 0.5
            match &= num_ids_considered - num_ids_matched <= 1
        else:
            match &= num_ids_matched == num_ids_considered
        # if not match:
        #     break
    match &= total_match_count>=1
    # if match:
    #     print(IDs)
    return match

# def medgenMatch(ClinVarIDs, infos):
#     count_medgen = 0 
#     count_medgen_correct = 0 
#     for cvID, info in zip(ClinVarIDs, infos):
#         if id.split(":")[0].lower() == "medgen" and usable_ClinVar_medgenID(info):
#             count_medgen+=1


def find_matching_indices(PIDs, disease):
    sdv_idx = [phenoInfo2bool(pids, disease) for pids in PIDs]
    return sdv_idx

# def find_matching_indices(PIDs, disease, Ontologies=Ontologies):
#     # v,c = np.unique(np.array([pid.split(':')[0] for pids in PIds for pid in re.split(r"[\,\;\|]", pids)]).flatten(), return_counts=True)
#     # print("PID counts:", dict(zip(v,c)))
#     #PIDs = [extractPIDs(s) for s in PIds]
#     Counts = {}
#     Counts_mapped_to_disease = {}
#     matches = {}
#     sdv_idx = [True for i in range(len(PIDs))]
#     for ont in Ontologies:
#         Counts[ont] = [len(pids) for pids in PIDs[ont]]
#         Counts_mapped_to_disease[ont] = [np.sum([find_in_disease_Phenotypes(pid, disease=disease) for pid in pids]) for pids in PIDs]  
#         if ont.lower() == "medgen":
#             matches[ont] = np.where(Counts[ont]<3, 
#                                     Counts[ont] - Counts_mapped_to_disease[ont]==0, 
#                                   Counts_mapped_to_disease[ont]/Counts[ont] > 0.5)
#         else:
#             matches[ont] = Counts[ont] - Counts_mapped_to_disease[ont]==0
#         sdv_idx &= matches[ont]

    # matches = [[phenoInfo2bool(info) for info in infos] for infos in Pheno_infos]
    #mdv_idx = [np.any(m) for m in matches]
    #sdv_idx = [np.all(m) if len(m) > 0 else False for m in matches]
#    return sdv_idx

def extractVariantsFromClinVar(clinvar_df, PIDs, disease):
    #mdv_idx, sdv_idx = find_matching_indices(PIDs, disease)
    sdv_idx = find_matching_indices(PIDs, disease)
    #MDVs = clinvar_df[mdv_idx]
    SDVs = clinvar_df[sdv_idx]
    #return MDVs, SDVs
    return SDVs

def concat(df_all, df, disease):
    print(f"Overlap {disease}:" + str(len(set(df.Variant.to_list()) & set(df_all.Variant.to_list()))))
    df_all = pd.concat([df_all, df], ignore_index=True)
    df_all1 = df_all.drop_duplicates()
    #assert len(df_all1) == len(df_all.drop_duplicates("Variant"))
    #df_all = df_all.groupby("Variant",as_index=False).agg({"GeneSymbol": "first", "Disease": lambda x: ",".join(x)})
    return df_all1
#def extractASDVariantsFromStudies(ASDVari, clinvar_df):

def extractRandomVariants(n, clinvar_df, exclude_variants):
    clinvar_df = clinvar_df[~clinvar_df.Variant.isin(exclude_variants)]
    P_idx, PLP_idx, _, BLB_idx, B_idx = clinSig_bool(clinvar_df.ClinicalSignificance)
    P_idx = P_idx | PLP_idx
    B_idx = B_idx | BLB_idx
    P_df = clinvar_df[P_idx]
    B_df = clinvar_df[B_idx]
    random_P = P_df.sample(n=n, random_state=42, replace=False)
    random_B = B_df.sample(n=n, random_state=42, replace=False)
    random_df = pd.concat([random_P, random_B], ignore_index=True)
    return random_df

# def filteredIDs(PIDs, Ontologies = Ontologies):
#     IDs = {ont:[] for ont in Ontologies}
#     for ids in PIDs:
#         for ont in Ontologies:
#             filtered_ids = []
#             for id in ids[ont]:
#                 info = find_in_ClinVar_master(id)
#                 if ont.lower() == "medgen":
#                     info = usable_ClinVar_medgenID(info)
#                 if info:
#                     filtered_ids.append(id)
#             IDs[ont].append(filtered_ids)
#     return IDs

def variantSummary(CV_df):
    nVars = CV_df.shape[0]
    nGenes = CV_df.GeneSymbol.nunique()
    P_idx, PLP_idx, VUS_idx, BLB_idx, B_idx = clinSig_bool(CV_df.ClinicalSignificance)
    print(f"Total variants: {nVars}")
    print(f"Total genes: {nGenes}")
    print(f"Pathogenic: {P_idx.sum()}")
    print(f"Pathogenic/Likely pathogenic: {PLP_idx.sum()}")
    print(f"Uncertain significance: {VUS_idx.sum()}")
    print(f"Benign: {BLB_idx.sum()}")
    print(f"Benign/Likely benign: {B_idx.sum()}")

if __name__ == "__main__":
    diseasesCV = ["CM","PCD", "PKD", "PH"]
    #diseasesCV = ["CM"]
    clinvar_df = pd.read_csv(ClinVar_file, sep='\t')
    clinvar_df = clinvar_df.rename(columns={"Variant_ID": "Variant"})
    PIDs = np.array([s for s in clinvar_df.PhenotypeIDS.to_list()])
    PIDs = [extractPhenoIDs(s)[1] for s in PIDs]
    #PIDs = filteredIDs(PIDs)
    # with open(f'PIDs.json', 'w') as f:
    #     json.dump(PIDs, f)
    # PIDs = onlyKeepDiseaseMedgenIDs(PIDs)
   
   
    for i, disease in enumerate(diseasesCV):
        print("\n\n")
        print(f'===================={disease}=====================')
        SDVs = extractVariantsFromClinVar(clinvar_df, PIDs, disease)
        #SDVs["Disease"] = d
        if i == 0:
            #MDVs_all = MDVs[["GeneSymbol","Variant_ID"]].copy()
            #SDVs_all = SDVs[["GeneSymbol","Variant_ID", "Disease"]].copy().reindex()
            SDVs_all = SDVs[["GeneSymbol","Variant"]].copy().reindex()
        else:
            #MDVs_all = pd.concat([MDVs_all, MDVs[["GeneSymbol","Variant_ID"]]], ignore_index=True).drop_duplicates()
            SDVs_all = concat(SDVs_all,SDVs[["GeneSymbol","Variant"]].copy().reindex(), disease)
        # print(f'#MDVs for {d}')
        # variantSummary(MDVs)        
        print(f'#SDVs for {disease}')
        variantSummary(SDVs)    
        filtered_folder = f'{DATA_DIR}/Filtered'
        #MDVs.to_csv(f'{filtered_folder}/{d}_MDVs.csv', index=False)
        with pd.ExcelWriter(extractedVariantsFile, mode="a", engine="openpyxl",  if_sheet_exists="replace") as writer:
            SDVs.to_excel(writer, sheet_name=f"{disease}", index=False)
        #SDVs.to_excel(extractedVariantsFile, sheet_name=f"{d}", index=False, mode="a", engine="openpyxl")

    ASDVars = pd.read_csv(ASDVarFile, sep='\t')
    ASDInCV = ASDVars.merge(clinvar_df, left_on=["GeneSymbol", "Variant"], right_on=["GeneSymbol", "Variant"], how='inner')[["GeneSymbol", "Variant", "ClinicalSignificance"]]
    ASDVars = ASDVars.merge(ASDInCV, left_on=["GeneSymbol", "Variant"], right_on=["GeneSymbol", "Variant"], how='left')
    ASDVars["In_ClinVar"] = ~ASDVars.ClinicalSignificance.isna()
    print(f"ASD variants in ClinVar: {ASDVars.In_ClinVar.sum()} out of {ASDVars.shape[0]}")
    print("Summary of ASD variants in ClinVar:")
    variantSummary(ASDInCV)
    with pd.ExcelWriter(extractedVariantsFile, mode="a", engine="openpyxl",  if_sheet_exists="replace") as writer:
        ASDVars.to_excel(writer, sheet_name=f"ASD", index=False)

    #SDVs_all = SDVs_all.rename(columns={"Variant_ID": "Variant"})
    SDVs_all = concat(SDVs_all, ASDVars[["GeneSymbol", "Variant"]].copy().reindex(), "ASD")
    SDVs_all["Target"] = True
    nonTarget_df = extractRandomVariants(10000, clinvar_df, SDVs_all.Variant.to_list())
    #nonTarget_df = nonTarget_df.rename(columns={"Variant_ID": "Variant"})
    with pd.ExcelWriter(extractedVariantsFile, mode="a", engine="openpyxl",  if_sheet_exists="replace") as writer:
        nonTarget_df.to_excel(writer, sheet_name="nonTarget", index=False)
    nonTarget_df["Target"] = False
    All_Vars = pd.concat([SDVs_all, nonTarget_df[["GeneSymbol", "Variant", "Target"]]], ignore_index=True)
    #assert All_Vars[["Variant", "GeneSymbol"]].duplicated().sum() == 0, "Duplicate variants found in All_Vars"
    with pd.ExcelWriter(extractedVariantsFile, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer:
        All_Vars.to_excel(writer, index=False, sheet_name="Variants_all")
    #All_Vars 

    #MDVs_all.to_csv(f'{filtered_folder}/All_MDV.csv', index=False)
   
from ExTSP.config import exTSP_file_gtex, exTSP_file_brainspan_cortical, ExtractedVariants_file, VariantSets, OLD_DATA_DIR
from ExTSP.commonFunctions import isCase_ASD, clinSig_bool
import pandas as pd
from pathlib import Path

def get_tripletSets_with_exTSP(disease, type="all", brainspan=False,  numVars=None):
    # type = "all", "P", "PLP", "B", "BLB", "VUS", "case", "control"
    # disease = "ASD", "CM", "PCD", "PKD", "PH", "nonTarget"
    if brainspan:
        exTSP_file = exTSP_file_brainspan_cortical
    else:
        exTSP_file = exTSP_file_gtex
    if not Path(exTSP_file).exists():
        exit("exTSP_file does not exist")
    df_extsp = pd.read_csv(exTSP_file, sep='\t')
    if not Path(ExtractedVariants_file).exists():
        exit("ExtractedVariants_file does not exist")
    df_vars = pd.read_excel(ExtractedVariants_file, sheet_name=disease)
    if type in ["P", "PLP", "B", "BLB", "VUS"]:
        P_idx, PLP_idx, VUS_idx, BLB_idx, B_idx = clinSig_bool(df_vars["ClinicalSignificance"])
        if type == "P":
            df_vars = df_vars[P_idx]
        elif type == "PLP":   
            df_vars = df_vars[PLP_idx]
        elif type == "B":
            df_vars = df_vars[B_idx]
        elif type == "BLB":   
            df_vars = df_vars[BLB_idx]
        elif type == "VUS":
            df_vars = df_vars[VUS_idx]  
        else:
            exit("Invalid type")
        df_vars = df_vars[["Variant", "ClinicalSignificance"]]
    elif type in ["case", "control"] and disease == "ASD":
        iscase = isCase_ASD(df_vars["Status"])
        if type == "case":
            df_vars = df_vars[iscase]
        elif type == "control":
            df_vars = df_vars[~iscase]
        else:
            exit("Invalid type")
        df_vars = df_vars[["Variant", "Status", "ClinicalSignificance", "In_ClinVar"]]
    df_vars.rename(columns={"ClinicalSignificance": "ClinVar_annotation"}, inplace=True)
    df_extsp = df_extsp.merge(df_vars, on="Variant", how="inner")
    if brainspan:
        df_extsp = df_extsp.drop_duplicates(subset=["Variant", "Transcript_id", "Developmental_Period", "Regioncode"])
    else:
        df_extsp = df_extsp.drop_duplicates(subset=["Variant", "Transcript_id", "Tissue"])
    if numVars is not None:
        df_extsp = sample_tripletSets(df_extsp, numVars)
    return df_extsp


def get_tripletSets_with_exTSP_old(disease, type="all", brainspan=False,  numVars=None):
    # type = "all", "P", "PLP", "B", "BLB", "VUS", "case", "control"
    # disease = "ASD", "CM", "PCD", "PKD", "PH", "nonTarget"
    diseases = ["CM", "PCD", "ASD"]
    nonTarget = False
    if disease not in diseases:
        disease = disease.split("_")[1]
        nonTarget = True
    exTSP_file_old = OLD_DATA_DIR / f"extsp_{disease}.tsv"
    exTSP_file_disease = OLD_DATA_DIR / f"classification/exTSP_file_expression_ensemble_annotation_{disease}.tsv"
    df_extsp_disease = pd.read_csv(exTSP_file_disease, sep='\t')
    if nonTarget:
        exTSP_file_non_disease = [OLD_DATA_DIR / f"classification/exTSP_file_expression_ensemble_annotation_{d}.tsv" for d in diseases if d != disease]
        df_extsp_non_disease = pd.DataFrame()
        for exTSP_file in exTSP_file_non_disease:
            if not Path(exTSP_file).exists():
                exit(f"exTSP_file {exTSP_file} does not exist")
            df_extsp_non_disease = pd.concat([df_extsp_non_disease, pd.read_csv(exTSP_file, sep='\t')])
        target_genes = df_extsp_disease.Gene.unique()
        df_extsp_non_disease = df_extsp_non_disease[~df_extsp_non_disease.Gene.isin(target_genes)]
        df_extsp = df_extsp_non_disease
    else:
        df_extsp = df_extsp_disease
    if type in ["P", "PLP", "B", "BLB", "VUS"] and disease != "ASD":
        P_idx, PLP_idx, VUS_idx, BLB_idx, B_idx = clinSig_bool(df_extsp["ClinVar_annotation"])
        if type == "P":
            df_extsp = df_extsp[P_idx]
        elif type == "PLP":   
            df_extsp = df_extsp[PLP_idx]
        elif type == "B":
            df_extsp = df_extsp[B_idx]
        elif type == "BLB":   
            df_extsp = df_extsp[BLB_idx]
        elif type == "VUS":
            df_extsp = df_extsp[VUS_idx]  
        else:
            exit("Invalid type")
    elif type in ["case", "control"] and disease == "ASD":
        iscase = isCase_ASD(df_extsp["Status"])
        if type == "case":
            df_extsp = df_extsp[iscase]
        elif type == "control":
            df_extsp = df_extsp[~iscase]
        else:
            exit("Invalid type")
    if numVars is not None:
        df_extsp = sample_tripletSets(df_extsp, numVars)
    return df_extsp



def sample_tripletSets(df_extsp, numVars):
    var_unique = df_extsp[["Variant"]].drop_duplicates()
    countRatio = int(numVars/len(var_unique))
    df_vars = pd.DataFrame({"Variant": []})
    if countRatio > 1:
        #df_vars = pd.DataFrame(var_unique*countRatio, columns=["Variant"])
        df_vars = pd.DataFrame(var_unique.Variant.unique().tolist()*countRatio, columns=["Variant"])
    df_vars2 = var_unique.sample(n=numVars%len(var_unique), replace=False)
    df_vars = pd.concat([df_vars, df_vars2], ignore_index=True)
    df_extsp_subsampled = df_vars.merge(df_extsp, on="Variant", how="inner")
    return df_extsp_subsampled

def variantCountsWithDuplicates(df):
    vCounts1 = df.Variant.value_counts()
    if "Tissue" in df.columns:
        vCounts2 = df[["Variant", "Transcript_id", "Tissue"]].drop_duplicates().Variant.value_counts()
    else:
        vCounts2 = df[["Variant", "Transcript_id", "Developmental_Period", "Regioncode"]].drop_duplicates().Variant.value_counts()
    vCounts3 = vCounts1/vCounts2
    return vCounts3.sum(), vCounts3


if __name__ == "__main__":
    for vset in VariantSets:
        if vset !="ASD":
            for type in ["P", "PLP", "B", "BLB", "VUS"]:
                df = get_tripletSets_with_exTSP(vset, type=type)
                print(f"=== {vset}_{type} ===")
                print(f"{df.Gene.nunique()} genes")
                print(f"{df.Variant.nunique()} variants")
                print(f"{len(df)} variant-transcript-tissue triplets")
        else:
            for type in ["case", "control"]:
                df = get_tripletSets_with_exTSP(vset, type=type)
                print(f"=== {vset}_{type} ===")
                print(f"{df.Gene.nunique()} genes")
                print(f"{df.Variant.nunique()} variants")
                print(f"{len(df)} variant-transcript-tissue triplets")
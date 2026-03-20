from ExTSP.config import exTSP_file, ExtractedVariants_file, VariantSets
from ExTSP.commonFunctions import isCase_ASD, clinSig_bool
import pandas as pd
from pathlib import Path

def get_tripletSets_with_exTSP(disease, type="all"):
    # type = "all", "P", "PLP", "B", "BLB", "VUS", "case", "control"
    # disease = "ASD", "CM", "PCD", "PKD", "PH", "nonTarget"
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
    df_extsp = df_extsp.drop_duplicates(subset=["Variant", "Transcript_id", "Tissue"])
    return df_extsp


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
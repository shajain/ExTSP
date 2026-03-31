from ExTSP.config import PTSE_file_gtex, PTSE_file_brainspan_cortical, MutPred2_processed_file
from ExTSP.config import exTSP_file_gtex, exTSP_file_brainspan_cortical, MANE_transcripts_file
import pandas as pd

from pathlib import Path

def generate_exTSP_scores(df_mp2, df_ptse, df_mane=None, brainspan=False):
    # Rename MANE transcript columns and create Transcript column containing ENST without version number
    if df_mane is not None:
        df_mane.rename(columns={"Ensembl_nuc": "Transcript_id"}, inplace=True)
        df_mane["Transcript"] = df_mane["Transcript_id"].str.split(".").str[0]
        df_mane = df_mane[["Transcript", "MANE_status"]]
        #Create a MANE_Select column with MANE_status == "MANE_Select" and MANE_Plus_Clinical column with MANE_status == "MANE_Plus_Clinical"
        df_mane["MANE_Select"] = df_mane["MANE_status"] == "MANE Select"
        df_mane["MANE_Plus_Clinical"] = df_mane["MANE_status"] == "MANE Plus Clinical"
        df_mane.drop(columns=["MANE_status"], inplace=True)
    # Select, Rename PTSE columns and create Transcript column containing ENST without version number
    if brainspan:
        df_ptse = df_ptse[["transcript_id","gene_name","Developmental_Period", "Regioncode","meanTPM","ptse","normalizedLR"]]
    else:
        df_ptse = df_ptse[["transcript_id","gene_name","tissue","meanTPM","ptse","normalizedLR"]]
    df_ptse.rename(columns={"transcript_id": "Transcript_id", "gene_name": "Gene", "tissue": "Tissue"}, inplace=True)
    df_ptse["Transcript"] = df_ptse["Transcript_id"].str.split(".").str[0]
    if df_mane is not None:
        # merge PTSE and MANE transcript files on Transcript column
        df_ptse = df_ptse.merge(df_mane, on="Transcript", how="left")
        df_ptse["MANE_Select"] = df_ptse["MANE_Select"].fillna(False)
        df_ptse["MANE_Plus_Clinical"] = df_ptse["MANE_Plus_Clinical"].fillna(False)
    # Select, Rename MutPred2 columns and create Transcript column containing ENST without version number
    df_mp2 = df_mp2[["Variant", "ENST",	"aa_sub","MutPred2","IsoPath"]]
    df_mp2["Transcript"] = df_mp2["ENST"].str.split(".").str[0]
    df_mp2.drop(columns=["ENST"], inplace=True)
    # Merge PTSE and MutPred2 files on Transcript column
    df_exTSP = df_mp2.merge(df_ptse, on=["Transcript"], how="inner")
    df_exTSP["exTSP"] = df_exTSP["IsoPath"] * df_exTSP["normalizedLR"]
    cols = ["Gene", "Transcript", "Transcript_id", "Variant", "MutPred2", "IsoPath", "ptse", "exTSP", "meanTPM"] 
    mane_cols = ["MANE_Select", "MANE_Plus_Clinical"]
    if brainspan:
        cols.extend(["Developmental_Period", "Regioncode"])
    else:
        cols.append("Tissue")
    if df_mane is not None:
        cols.extend(mane_cols)
    df_exTSP = df_exTSP[cols]
    return df_exTSP

if __name__ == "__main__":
    if not Path(PTSE_file_gtex).exists():
        exit("Run compute_PTSE.py first")
    if not Path(MutPred2_processed_file).exists():
        exit("Run mutPred2IsoPath.py first")
    # Read PTSE, MutPred2 processed files, and MANE transcript files
    if not Path(exTSP_file_gtex).exists():
        df_ptse = pd.read_csv(PTSE_file_gtex, sep='\t')
        df_mp2 = pd.read_csv(MutPred2_processed_file, sep='\t')
        df_mane = pd.read_csv(MANE_transcripts_file, sep = "\t")
        df_exTSP = generate_exTSP_scores(df_mp2, df_ptse, df_mane)
        df_exTSP.to_csv(exTSP_file_gtex, sep='\t', index=False)
        print(f"Saved exTSP file to {exTSP_file_gtex}")
    if not Path(PTSE_file_brainspan_cortical).exists():
        exit("Run compute_PTSE.py first")    
    if not Path(exTSP_file_brainspan_cortical).exists():
        df_ptse = pd.read_csv(PTSE_file_brainspan_cortical, sep='\t')
        df_ptse.rename(columns={"transcript": "transcript_id", "external_gene_name": "gene_name",}, inplace=True)
        df_mp2 = pd.read_csv(MutPred2_processed_file, sep='\t')
        df_mane = pd.read_csv(MANE_transcripts_file, sep = "\t")
        df_exTSP = generate_exTSP_scores(df_mp2, df_ptse, df_mane, brainspan=True)
        df_exTSP.to_csv(exTSP_file_brainspan_cortical, sep='\t', index=False)
        print(f"Saved exTSP file to {exTSP_file_brainspan_cortical}")

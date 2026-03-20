from ExTSP.config import PTSE_file, MutPred2_processed_file, exTSP_file, MANE_transcripts_file
import pandas as pd

from pathlib import Path



if __name__ == "__main__":
    if not Path(PTSE_file).exists():
        exit("Run GTEx2PTSE.py first")
    if not Path(MutPred2_processed_file).exists():
        exit("Run mutPred2IsoPath.py first")
    # Read PTSE, MutPred2 processed files, and MANE transcript files
    df_ptse = pd.read_csv(PTSE_file, sep='\t')
    df_mp2 = pd.read_csv(MutPred2_processed_file, sep='\t')
    mane_transcript = pd.read_csv(MANE_transcripts_file, sep = "\t")
    # Rename MANE transcript columns and create Transcript column containing ENST without version number
    mane_transcript.rename(columns={"Ensembl_nuc": "Transcript_id"}, inplace=True)
    mane_transcript["Transcript"] = mane_transcript["Transcript_id"].str.split(".").str[0]
    mane_transcript = mane_transcript[["Transcript", "MANE_status"]]
    #Create a MANE_Select column with MANE_status == "MANE_Select" and MANE_Plus_Clinical column with MANE_status == "MANE_Plus_Clinical"
    mane_transcript["MANE_Select"] = mane_transcript["MANE_status"] == "MANE Select"
    mane_transcript["MANE_Plus_Clinical"] = mane_transcript["MANE_status"] == "MANE Plus Clinical"
    mane_transcript.drop(columns=["MANE_status"], inplace=True)
    # Select, Rename PTSE columns and create Transcript column containing ENST without version number
    df_ptse = df_ptse[["transcript_id","gene_name","tissue","meanTPM","ptse","normalizedLR"]]
    df_ptse.rename(columns={"transcript_id": "Transcript_id", "gene_name": "Gene",
                    "tissue": "Tissue"}, inplace=True)
    df_ptse["Transcript"] = df_ptse["Transcript_id"].str.split(".").str[0]
    # merge PTSE and MANE transcript files on Transcript column
    df_ptse = df_ptse.merge(mane_transcript, on="Transcript", how="left")
    df_ptse["MANE_Select"] = df_ptse["MANE_Select"].fillna(False)
    df_ptse["MANE_Plus_Clinical"] = df_ptse["MANE_Plus_Clinical"].fillna(False)
    # Select, Rename MutPred2 columns and create Transcript column containing ENST without version number
    df_mp2 = df_mp2[["Variant", "ENST",	"aa_sub","MutPred2","IsoPath"]]
    df_mp2["Transcript"] = df_mp2["ENST"].str.split(".").str[0]
    df_mp2.drop(columns=["ENST"], inplace=True)
    # Merge PTSE and MutPred2 files on Transcript column
    df_exTSP = df_mp2.merge(df_ptse, on=["Transcript"], how="left")
    df_exTSP["exTSP"] = df_exTSP["IsoPath"] * df_exTSP["normalizedLR"]
    df_exTSP = df_exTSP[["Gene", "Transcript", "Transcript_id", "Variant", "Tissue",  
                        "MutPred2", "IsoPath", "ptse", "exTSP", "meanTPM", "MANE_Select", "MANE_Plus_Clinical"]]
    df_exTSP.to_csv(exTSP_file, sep='\t', index=False)
    print(f"Saved exTSP file to {exTSP_file}")

   # "Gene	Transcript	Transcript_id	Variant	ClinVar_annotation	Tissue	IsoPath	MutPred2	ptse	exTSP	MANE_Select	MANE_Plus_Clinical	meanTPM"



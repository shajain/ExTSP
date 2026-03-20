from ExTSP.config import MutPred2_filtered_file, MutPred2PosteriorTable_file, MutPred2_processed_file
from ExTSP.config import VariantAnnotations_file, MANE_transcripts_file
import pandas as pd
#from Bio import SeqIO
from pathlib import Path
import numpy as np
def interpolate_from_table(values: np.ndarray, table: pd.DataFrame) -> np.ndarray:
    """
    Linearly interpolates values from a two-column lookup table.

    Args:
        values:  1D array of query values to look up
        table:   DataFrame with two or more columns — first column is the x (key), second column is the y (value)

    Returns:
        1D array of interpolated y values
    """
    x_col, y_col = table.columns[0], table.columns[1]
    table_sorted = table.sort_values(x_col).drop_duplicates(subset=x_col)

    return np.interp(
        values,
        table_sorted[x_col].values,
        table_sorted[y_col].values
    )

def main():
    if not Path(MutPred2_filtered_file).exists():
        exit("MutPred2_filtered_file does not exist")
    df_mp2 = pd.read_csv(MutPred2_filtered_file, sep='\t')
    df_mp2 = df_mp2[~df_mp2["Substitution"].isna()]
    if not Path(MutPred2PosteriorTable_file).exists():
        exit("MutPred2PosteriorTable_file does not exist")
    df_mp2_posterior = pd.read_csv(MutPred2PosteriorTable_file, sep='\t')
    if not Path(VariantAnnotations_file).exists():
        exit("VariantAnnotations_file does not exist")
    df_VarAnnotations = pd.read_csv(VariantAnnotations_file, sep='\t')
    df_VarAnnotations = df_VarAnnotations.drop_duplicates()
    df_VarAnnotations.drop(columns=["chr", "pos", "ref", "alt"], inplace=True)
    print(f"Number of variant-transcript pairs in VariantAnnotations: {len(df_VarAnnotations)}")
    df_mp2_merged = df_mp2.merge(df_VarAnnotations, left_on=["ID","Substitution"], right_on=["ENST","aa_sub"], how="left")
    df_mp2_merged.drop(columns=["ID","Substitution"], inplace=True)
    df_mp2_merged = df_mp2_merged[["Variant", "ENST", "ENSG", "ENSP", "gene", "aa_sub", "MutPred2"]]
    assert not df_mp2_merged.Variant.isna().any(), "Some MutPred2 entires not mapped back to Genomic variants"
    # df_mp2_IsoPath = df_mp2_merged.merge(df_mp2_posterior, right_on=["score"], left_on=["Thresholds"], how="left")
    df_mp2_merged["IsoPath"] = interpolate_from_table(df_mp2_merged["MutPred2"].values, df_mp2_posterior)
    print(f"Number of variants with MP2 Scores: {df_mp2_merged.Variant.nunique()}")
    print(f"Number of variant-transcript pairs with MP2 Scores: {len(df_mp2_merged)}")
    df_mp2_merged.to_csv(MutPred2_processed_file, sep='\t', index=False)
    print(f"Saved MutPred2 processed file to {MutPred2_processed_file}")


if __name__ == "__main__":
    main()
    
    
        
    

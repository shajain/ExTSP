from ExTSP.config import MutPred2_gtex_file, MutPred2PosteriorTable_file, IsoPath_gtex_file,MutPred2_brainspan_file, IsoPath_brainspan_file
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

def add_isoPath_column(df_mp2):
    if not Path(MutPred2PosteriorTable_file).exists():
        exit("MutPred2PosteriorTable_file does not exist")
    df_mp2_posterior = pd.read_csv(MutPred2PosteriorTable_file, sep='\t')
    df_mp2["IsoPath"] = interpolate_from_table(df_mp2["MutPred2"].values, df_mp2_posterior)
    return df_mp2

def main(mutPred2_file, isoPath_file):
    if not Path(mutPred2_file).exists():
        exit("mutPred2_file does not exist")
    df_mp2 = pd.read_csv(mutPred2_file, sep='\t')
    df_mp2 = df_mp2[~df_mp2["aa_sub"].isna()]
    df_mp2 = df_mp2[~df_mp2["MutPred2"].isna()]
    assert not df_mp2["Variant"].isna().any(), "Some MutPred2 entires not mapped back to Genomic variants"
    # if "Variant" in df_mp2.columns:
    #     df_mp2.rename(columns={"Variant": "Variant_old"}, inplace=True)

    # if not Path(VariantAnnotations_file).exists():
    #     exit("VariantAnnotations_file does not exist")
    # df_VarAnnotations = pd.read_csv(VariantAnnotations_file, sep='\t')
    # df_VarAnnotations = df_VarAnnotations.drop_duplicates()
    # df_VarAnnotations.drop(columns=["chr", "pos", "ref", "alt"], inplace=True)
    # print(f"Number of variant-transcript pairs in VariantAnnotations: {len(df_VarAnnotations)}")
    # df_mp2_merged = df_mp2.merge(df_VarAnnotations, left_on=["ID","Substitution"], right_on=["ENST","aa_sub"], how="left")

    # if "Variant_old" in df_mp2_merged.columns:
    #     #is_not_na_variant = np.logical_and(df_mp2_merged["Variant_old"].notna(), df_mp2_merged["Variant"].notna())
    #     #assert df_mp2_merged.loc[is_not_na_variant].Variant.equals(df_mp2_merged.loc[is_not_na_variant].Variant_old), "Variant and Variant_old are not the same"
    #     is_na_variant = df_mp2_merged["Variant"].isna()
    #     df_mp2_merged.loc[is_na_variant, "Variant"] = df_mp2_merged.loc[is_na_variant, "Variant_old"]
    #     df_mp2_merged.drop(columns=["Variant_old"], inplace=True)
    #     print(f"Number of variants not mapped back to Genomic variants: {len(df_mp2_merged[df_mp2_merged['Variant'].isna()])}")

    # df_mp2_merged.drop(columns=["ID","Substitution"], inplace=True)
    # df_mp2_merged = df_mp2_merged[["Variant", "ENST", "ENSG", "ENSP", "gene", "aa_sub", "MutPred2"]]
    #assert not df_mp2_merged.Variant.isna().any(), "Some MutPred2 entires not mapped back to Genomic variants"
    # df_mp2_IsoPath = df_mp2_merged.merge(df_mp2_posterior, right_on=["score"], left_on=["Thresholds"], how="left")
    # df_mp2_merged["IsoPath"] = interpolate_from_table(df_mp2_merged["MutPred2"].values, df_mp2_posterior)
    #df_mp2_merged = add_isoPath_column(df_mp2_merged)
    df_mp2_isoPath = add_isoPath_column(df_mp2)
    print(f"Number of variants with MP2 Scores: {df_mp2_isoPath.Variant.nunique()}")
    print(f"Number of variant-transcript pairs with MP2 Scores: {len(df_mp2_isoPath)}")
    df_mp2_isoPath.to_csv(isoPath_file, sep='\t', index=False)
    print(f"Saved MutPred2 processed file to {isoPath_file}")


if __name__ == "__main__":
    #main(mutPred2_file=MutPred2_gtex_file, isoPath_file=IsoPath_gtex_file)
    main(mutPred2_file=MutPred2_brainspan_file, isoPath_file=IsoPath_brainspan_file)
    
    
        
    

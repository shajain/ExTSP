from ExTSP.config import VEP_output_DIR, VariantAnnotations_file, ExtractedVariants_file
import pandas as pd
# from Bio import SeqIO
from pathlib import Path


if __name__ == "__main__":
    #Its an excel file
    df_vars = pd.read_excel(ExtractedVariants_file, sheet_name="Variants_all")
    for chr in range(1,25):
        df = pd.read_csv(VEP_output_DIR / f"dbNSFP5.3a_variant.chr{chr}.filtered.vep.txt", sep='\t', dtype=str)
        df["Variant"] = df["chr"] + "_" + df["pos"] + "_" + df["ref"] + "/" + df["alt"]
        df_filtered = df[df["Variant"].isin(df_vars["Variant"])]
        if chr==1:
            DF = df_filtered
        else:
            DF = pd.concat([DF, df_filtered])
        print(f"read chr{chr}, Number of mapped variant-transcripts: {len(DF)}")
    print(f"Number of Extracted Variants: {df_vars.Variant.nunique()}")
    print(f"Number of Extracted Variants mapped to sinlge AA substitutions: {DF.Variant.nunique()}")
    DF.to_csv(VariantAnnotations_file, sep='\t', index=False)
    print(df.head())

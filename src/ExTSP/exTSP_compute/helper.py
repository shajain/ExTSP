from ExTSP.config import VariantAnnotations_file
import pandas as pd

def filter_by_transcripts(df,  transcripts, transcript_column = "transcript_id",):
    transcripts = [tr.split('.')[0] for tr in transcripts]
    df["transcript"] = df[transcript_column].str.split('.').str[0]
    df = df[df["transcript"].isin(transcripts)]
    return df


def map_to_genomic_variants(df, transcript_column = "ID", aa_substitution_column = "Substitution"):
    df_annotations = pd.read_csv(VariantAnnotations_file, sep="\t")
    df["transcript"] = df[transcript_column].str.split('.').str[0]
    df_annotations["transcript"] = df_annotations["ENST"].str.split('.').str[0]
    df1 = df.merge(df_annotations, left_on=["transcript", aa_substitution_column], right_on=["transcript", "aa_sub"], how="left")
    pairs_in_df = set(df[["transcript", aa_substitution_column]].to_records(index=False).tolist())
    pairs_in_df1 = set(df1[["transcript", "aa_sub"]].to_records(index=False).tolist())
    not_in_annotations = len(pairs_in_df - pairs_in_df1)
    print(f"number of variants not in annotations: {not_in_annotations}")
    return df1
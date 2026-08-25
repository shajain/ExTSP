# %% IMPORT STATEMENTS
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP
from ExTSP.commonFunctions import clinSig_bool, isCase_ASD

# %% LOAD TIM'S MATCHING
import pandas as pd
data_path = "/Users/karnamendonca-kamath/Desktop/exTSP2/gtex_brainspan_transcripts_mapped_to_ensembl_v113.txt"
mapping = pd.read_csv(data_path, sep='\t')

### 387 rows have nulls in the sequence columns
### NO other nulls anywhere
mapping

# %%
mapping['transcript_v113'].nunique()


# %%
### PART 2: For each disease, load the data.

def summarize_sets(df, disease):
    """Collect the unique Genes/Transcripts/Variants (overall and by ClinVar/case-control category) in df."""
    genes = set(df['Gene'].dropna().unique())
    transcripts = set(df['Transcript_id'].dropna().unique())

    if disease == "ASD":
        variants_df = df[['Variant', 'Status']].drop_duplicates(subset='Variant')
        variants = set(variants_df['Variant'])

        iscase = isCase_ASD(variants_df['Status'])
        case_variants = set(variants_df.loc[iscase, 'Variant'])
        control_variants = set(variants_df.loc[~iscase, 'Variant'])

        return {
            "Number of Genes": genes,
            "Number of Total Transcripts": transcripts,
            "Number of Total Variants": variants,
            "Number of Case Variants": case_variants,
            "Number of Control Variants": control_variants,
        }

    variants_df = df[['Variant', 'ClinVar_annotation']].drop_duplicates(subset='Variant')
    variants = set(variants_df['Variant'])

    _, PLP_idx, VUS_idx, BLB_idx, _ = clinSig_bool(variants_df['ClinVar_annotation'])
    plp_variants = set(variants_df.loc[PLP_idx, 'Variant'])
    blb_variants = set(variants_df.loc[BLB_idx, 'Variant'])
    vus_variants = set(variants_df.loc[VUS_idx, 'Variant'])

    return {
        "Number of Genes": genes,
        "Number of Total Transcripts": transcripts,
        "Number of Total Variants": variants,
        "Number of P/LP Variants": plp_variants,
        "Number of B/LB Variants": blb_variants,
        "Number of VUS": vus_variants,
    }


def print_stats(stats, baseline=None):
    """Print each stat's count; if baseline is given, also print how many/what % were lost relative to it."""
    for label, current_set in stats.items():
        n = len(current_set)
        if baseline is None:
            print(f"  {label}: {n}")
        else:
            base_set = baseline[label]
            n_lost = len(base_set - current_set)
            pct_lost = 100 * n_lost / len(base_set) if base_set else 0.0
            print(f"  {label}: {n} ({n_lost} lost, {pct_lost:.2f}% of {len(base_set)})")


for disease in ["CM", "PCD", "ASD", "PKD", "PH", "nonTarget"]:
    print(f"=== {disease} ===")
    df = get_tripletSets_with_exTSP(disease)

    ### First print the dataset stats
    # Number of Genes
    # Number of Total Transcripts
    # Number of Total Variants
    # Number of P/LP Variants
    # Number of B/LB Variants
    # Number of VUS
    stats_before = summarize_sets(df, disease)
    print("Before ensembl v113 mapping:")
    print_stats(stats_before)

    current_mapping = mapping[mapping['dataset'] == 'gtex']

    ### Then left join this with the current_mapping dataframe on transcript_id
    # df's Transcript_id is versioned (e.g. ENST00000359218.10) for the gtex exTSP file,
    # which matches current_mapping's versioned 'transcript_id' column (not the
    # unversioned 'transcript' column used for the BrainSpan join above).
    merged = df.merge(
        current_mapping,
        left_on="Transcript",
        right_on="transcript",
        how="left",
        suffixes=("", "_mapping"),
    )
    matched = merged[merged["transcript"].notna()]

    # Print all the same stats, and the percentage and number of variants lost in parenthesis
    stats_after = summarize_sets(matched, disease)
    print("After left join with v113 mapping (gtex):")
    print_stats(stats_after, baseline=stats_before)
    print()
# %%

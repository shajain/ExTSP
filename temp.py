# %%
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP
import pandas as pd
from ExTSP.classification.auc_bootstrap import remove_conflicting_variants_v2
pd.set_option('display.max_columns', 500)

# %%
conflicting = pd.read_csv("/home/karnapmk/ExTSP/src/ExTSP/data/summary_of_conflicting_interpretations.txt", sep="\t")
conflicting
# %%
for disease in ['PCD', 'CM', 'nonTarget']:
    for subset in ['P', 'PLP', 'B', 'BLB', 'VUS']:
        df = get_tripletSets_with_exTSP(disease, subset)
        print(f"Number of {disease} {subset} Variants before filtering: {len(df.drop_duplicates(['Gene', 'Variant']))}")
        df = remove_conflicting_variants_v2(df, conflicting)
        print(f"Number of {disease} {subset} Variants after filtering: {len(df.drop_duplicates(['Gene', 'Variant']))}")
# %%

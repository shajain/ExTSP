Run from the root folder
pip install -e .

To get the datasets for a disease Run
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP
df_target_all = get_tripletSets_with_exTSP("CM")
df_target_plp = get_tripletSets_with_exTSP("CM", "PLP")
df_target_blb = get_tripletSets_with_exTSP("CM", "BLB")
df_nonTarget_plp = get_tripletSets_with_exTSP("nonTarget", "PLP")
df_nonTarget_blb = get_tripletSets_with_exTSP("nonTarget", "BLB")


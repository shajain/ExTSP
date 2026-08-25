import sys
sys.path.insert(0, "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/scripts")
from fig_3b import load_vus, select_gene_and_variants

locked_genes = {
    "CM": "TPM1",
    "PCD": "DNAAF3",
    "ASD": "SCN2A",
}
annotations = {
    "CM": "Uncertain significance",
    "PCD": "Uncertain significance",
    "ASD": "ASD",
}

for disease, gene in locked_genes.items():
    sub_df = load_vus(disease, annotations[disease])
    try:
        g, vt_pairs, hit_variants = select_gene_and_variants(sub_df, gene_pool=[gene])
    except ValueError as e:
        print(f"{disease} ({gene}): {e}")
        continue
    print(f"\n{disease}: gene={g}")
    print(f"  vt_pairs = {vt_pairs}")
    print(f"  hit_variants = {hit_variants}")

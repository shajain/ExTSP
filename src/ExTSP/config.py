from importlib.resources import files


NUM_CPU = 7
NBOOT = 100

#DATA_DIR = "data"
DATA_DIR = files("ExTSP") / "data"
RESULTS_DIR = files("ExTSP") / "results"
FIGURES_DIR = files("ExTSP") / "figures"
RAW_DIR = DATA_DIR / "raw"
FILTERED_DIR = DATA_DIR / "filtered"
PROCESSED_DIR = DATA_DIR / "processed"
PHENOTYPES_DIR = DATA_DIR / "phenotypes"
ClinVar_master_file = RAW_DIR / "missense_variants_ClinVar_Nov_2024.tsv"
ClinVar_phenotypes_master_file = PHENOTYPES_DIR / "ClinVar_Phenotypes_Master.json"
ExtractedVariants_file = FILTERED_DIR / "extractedVariants.xlsx"
MutPred2_filtered_file = FILTERED_DIR / "MutPred2_for_filtered_variants.tsv"
MutPred2_processed_file = PROCESSED_DIR / "MutPred2_processed_with_IsoPath.tsv"
MutPred2PosteriorTable_file = RAW_DIR / "MutPred2_local_posteriors_new.tsv"
MANE_transcripts_file = RAW_DIR / "MANE.GRCh38.v1.3.summary.tsv"
GTeX_file = RAW_DIR / "GTeX_V8_short_read_transcripts_meanTPM_expression_across_tissue_samples_Gencodev26_without_tpm_filter.tsv"
BrainSpan_expression_file_cortical = RAW_DIR / "BrainSpan_FuSatt_meanTPM_byRegion_byPeriod_Cortical_Chau11.tsv"
PTSE_file_gtex = PROCESSED_DIR / "GTEx_PTSE.tsv"
PTSE_file_brainspan_cortical = PROCESSED_DIR / "BrainSpan_PTSE_cortical.tsv"
Fasta_DIR = FILTERED_DIR / "fastachromosomes"
VEP_output_DIR = FILTERED_DIR / "gtex_veps_filtered"
VariantAnnotations_file = FILTERED_DIR / "variant_annotations.tsv"
exTSP_file_gtex = PROCESSED_DIR / "exTSP_gtex.tsv"
exTSP_file_brainspan_cortical = PROCESSED_DIR / "exTSP_brainspan_cortical.tsv"




DISEASE_SEMANTIC_TYPES = ["disease or syndrome", "congenital abnormality", "neoplastic process", "mental or behavioral dysfunction", "pathologic function",
                              "anatomical abnormality", "finding"]

NCBI_API_KEY = "1b42f5985b629bc2dd1cb5016ccae752db09"

OMIM_API_KEY = "XW0cWbGfTiaaarn0vEAASw"

Diseases = {"CM": "Cardiomyopathy", 
            "PKD": "Polycystic Kidney Disease", 
            "PH": "Pulmonary Hypertension", 
            "PCD": "Primary Ciliary Dyskinesia", 
            "ASD": "Autism Spectrum Disorder"}

VariantSets = list(Diseases.keys()) + ["nonTarget"]

exTSP_disease_files = {}
for disease in Diseases:
    exTSP_disease_files[disease] = f"extTSP_{disease}.tsv"
exTSP_disease_files["nonTarget"] = "extTSP_Target.tsv"
#Diseases = {"CM": "Cardiomyopathy"}

Ontologies = ["OMIM", "MONDO", "MedGen"]

Phenotype_search_parameters_OMIM = {"CM": {"search_str": "cardiomyopathy"},
                                    "PCD": {"search_str": "primary ciliary dyskinesia"},    
                                    "PKD": {"search_str": "polycystic kidney disease"},  
                                    "PH": {"search_str": "pulmonary hypertension"}}

Phenotype_search_parameters_MONDO = {"CM": {"search_str": "MONDO:0004994"},
                                    "PCD": {"search_str": "MONDO:0016575"},
                                    "PKD": {"search_str": "MONDO:0020642"},
                                    "PH": {"search_str": "MONDO:0005149"}}


Phenotype_search_parameters_MedGen = {"CM": {"search_str": "cardiomyopathy", "filter_terms": ["cardiomyopathy", "ventricular", "myocardial", 
                                            "arrhythmia", "heart failure", "cardiomegaly", "cardiac"]},
                                    "PCD": {"search_str": "primary ciliary dyskinesia", "filter_terms": ["ciliary", "cilia"]},
                                    "PKD": {"search_str": "polycystic kidney disease", "filter_terms": ["kidney", "renal", "polycystic" ,"nephro"]},
                                    "PH": {"search_str": "pulmonary hypertension", "filter_terms": ["pulmonary", "hypertension", "lung", 
                                                                                "arterial", "pulmonic", "right heart"]}}
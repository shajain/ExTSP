import pronto
from pathlib import Path
from ExTSP.config import DATA_DIR
from ExTSP.Phenotypes.common_functions import find_in_disease_Phenotypes, find_in_ClinVar_master

Phenotypes_folder = f"{DATA_DIR}/phenotypes"
clinvar_phenotypes_file = f"{Phenotypes_folder}/clinvar_MONDO_master.json"

ONT = pronto.Ontology("data/raw/mondo.obo")


def find_disease_Phenotypes_MONDO(search_params):   
    mondo_id = search_params["search_str"]
    if not mondo_id.startswith("MONDO"):    
        mondo_id = f"MONDO:{mondo_id}"
    term = ONT[mondo_id]  
    descendants = list(term.subclasses(distance=None))
    IDs = {d.id:{"title": d.name, "alternative_titles":  "|".join([s.description for s in d.synonyms])} for d in descendants}
    alt_titles = "|".join([s.description for s in term.synonyms])
    IDs[term.id] = {"title": term.name, "alternative_titles": alt_titles}  # include the original term
    # if mondo_id not in IDs:
    #     print(f"Warning: Mondo ID {mondo_id} not found in ontology. Adding it with the provided disease name.")
    #     IDs[mondo_id] = {"title": disease, "alternative_titles": ""}      # include the original term
    return IDs

def MONDO_extraInfo(mondo_id):
    if not mondo_id.startswith("MONDO"):    
        mondo_id = f"MONDO:{mondo_id}"
    if mondo_id in ONT:
        term = ONT[mondo_id]  
        alt_titles = "|".join([s.description for s in term.synonyms])
        return {"title": term.name, "alternative_titles": alt_titles}   
    return None
    
    

def get_clinVar_MONDO_master(mondo_ids):
    mondo_ids_info = {id: MONDO_extraInfo(id) for id in mondo_ids}
    mondo_ids_info_filtered = {id: info for id, info in mondo_ids_info.items() if info}
    # with open(clinvar_phenotypes_file, "w") as f:
    #         json.dump(mondo_ids_info_filtered, f)
    #         print(f"Saved {len(mondo_ids_info_filtered)} filtered MONDO IDs to file")
    return mondo_ids_info_filtered

# def get_clinVar_MONDO_master(mondo_ids):
#     file_path = Path(clinvar_phenotypes_file)
#     if file_path.exists():
#         with open(clinvar_phenotypes_file, "r") as f:
#             mondo_ids_info_filtered = json.load(f)
#             print(f"Loaded {len(mondo_ids_info_filtered)} filtered MONDO IDs from file")
#     else:
#         mondo_ids_info_filtered = create_clinVar_MONDO_master(mondo_ids)
#     return mondo_ids_info_filtered


if __name__ == "__main__":
    Ids = find_disease_Phenotypes_MONDO({"search_str": "MONDO:0004994"})
    print(Ids)
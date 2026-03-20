from ExTSP.config import DATA_DIR, Ontologies
import json
from pathlib import Path
import requests
import time


disease_phenotype_file = f"{DATA_DIR}/Phenotypes/Disease_Phenotypes.json"
clinvar_master_file = f"{DATA_DIR}/Phenotypes/ClinVar_Phenotypes_Master.json"
file_path = Path(disease_phenotype_file)
PHENOTYPES_SEARCH = {}
CLINVAR_MASTER = {}

def find_in_disease_Phenotypes(id, disease=None, ontology=None):
    ontologies = [ontology] if ontology else Ontologies
    split = id.split(":")
    if len(split)==1:
        if ontology:
            id = f"{ontology}:{id}"
        else:
            raise ValueError(f"Warning: ID {id} does not have an ontology prefix. Consider adding one for more accurate searching.")
    elif (len(split) == 2 and split[0] not in ontologies) or len(split)>2:
        print(f"Not a valid ontology or ontology not supported:{id}")        

    PHENOTYPES_SEARCH = load_disease_Phenotypes()
    if PHENOTYPES_SEARCH:
        diseases = [disease] if disease else list(PHENOTYPES_SEARCH.keys())
        phenotypes_search_results = {}
        for d in diseases:
            for o in ontologies:
                phenotypes_search_results |= PHENOTYPES_SEARCH[d].get(o, {}).get("search_results", {}) 
        if id in phenotypes_search_results:
            return phenotypes_search_results[id]
    return {}


def find_in_ClinVar_master(id, ontology=None):
    if len(id.split(":"))==1:
        if ontology:
            id = f"{ontology}:{id}"
        else:
            raise ValueError(f"Warning: ID {id} does not have an ontology prefix. Consider adding one for more accurate searching.")
    CLINVAR_MASTER = load_clinVar_master()
    if id in CLINVAR_MASTER:
        return CLINVAR_MASTER[id]
    return {}


def load_disease_Phenotypes():
    global PHENOTYPES_SEARCH
    if not PHENOTYPES_SEARCH:
        file_path = Path(disease_phenotype_file)
        if file_path.exists():
            with open(file_path, "r") as f:
                PHENOTYPES_SEARCH = json.load(f)
    return PHENOTYPES_SEARCH

def load_clinVar_master():
    global CLINVAR_MASTER
    if not CLINVAR_MASTER:
        file_path = Path(clinvar_master_file)
        if file_path.exists():
            with open(file_path, "r") as f:
                CLINVAR_MASTER = json.load(f)
    return CLINVAR_MASTER


def get_response(url, params):
    success = False
    while not success:
        response = requests.get(url, params=params)
        try:
            response.raise_for_status()
            success = True
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:  # Too Many Requests
            #print("Rate limit hit. Sleeping...")
                time.sleep(0.5)  # Sleep for a short time before retrying
    return response
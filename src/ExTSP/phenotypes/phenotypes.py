#from urllib import response
import json
import pandas as pd
from pathlib import Path
from ExTSP.config import RAW_DIR, PHENOTYPES_DIR, ClinVar_master_file, ClinVar_phenotypes_master_file
from ExTSP.config import DATA_DIR, Diseases, Ontologies
from ExTSP.config import Phenotype_search_parameters_MedGen, Phenotype_search_parameters_OMIM, Phenotype_search_parameters_MONDO
from ExTSP.Phenotypes.omim import get_clinVar_OMIM_master, find_disease_Phenotypes_OMIM
from ExTSP.Phenotypes.mondo import get_clinVar_MONDO_master, find_disease_Phenotypes_MONDO
from ExTSP.Phenotypes.medgen import find_disease_Phenotypes_MedGen, get_clinVar_MedGen_master
from ExTSP.extractVaraints.usefulFuncs import extractPhenoIDs





# def search_omim(disease, start=0, limit=100):
#     API_KEY = "XW0cWbGfTiaaarn0vEAASw"
#     BASE_URL = "https://api.omim.org/api"

#     params = {
#         "search": f"title:{disease}",
#         "format": "json",
#         "apiKey": API_KEY,
#         "start": start,
#         "limit": limit,  
#     }
#     IDs = {}
#     while True:
#         params["start"] = start
#         response = requests.get(f"{BASE_URL}/entry/search", params=params)
#         data = response.json()
    
#         if response.status_code == 200:
#             response = response.json()
#         else:
#             break
#         entries = data['omim']['searchResponse']['entryList']
#         for entry in entries:
#             e = entry['entry']
#             #print(e['mimNumber'], e['titles']['preferredTitle'])
#             IDs['OMIM:'+str(e['mimNumber'])] = {"title": e['titles']['preferredTitle']}
#         start += len(entries)
#         if start >= data['omim']['searchResponse']['totalResults']:
#             break
#     return IDs

# def search_hpo(hpo_id):
#     url = f"https://ontology.jax.org/api/hp/terms/{hpo_id}/descendants"
#     response = requests.get(url)

#     if response.status_code == 200:
#         data = response.json()
#     else:
#         print(f"Error: {response.status_code}")
#         return []
#     IDs = [hpo_id]
#     for entry in data:
#         id = entry['id']
#         IDs.append(id)
#     return IDs






def getDiseaseIDs(disease):
    with open(f"{PHENOTYPES_DIR}/{disease}_IDs.json", "r") as f:
        data = json.load(f)
    medgen_ids = {k:v for k,v in data["medgen_ids"].items() if v["hierarchy"] != "ancestor"}
    IDs = set(list(data["omim_ids"].keys())  + list(data["mondo_ids"].keys()) + list(medgen_ids.keys()))
    return IDs
    
if __name__ == "__main__":
    Phenotype_search_results = {}
    for disease in Diseases.keys():
        if disease == "ASD":
            continue
        Phenotype_search_results[disease] = {}
        Phenotype_search_results[disease]["OMIM"] = {}
        Phenotype_search_results[disease]["MONDO"] = {}
        Phenotype_search_results[disease]["MedGen"] = {}
        file = f'{PHENOTYPES_DIR}/{disease}_Phenotypes.json'
        file_path = Path(file)
        # if file_path.exists():
        #     continue
        print(f"Searching for {disease} in OMIM...")
        omim_search_results = find_disease_Phenotypes_OMIM(Phenotype_search_parameters_OMIM[disease])
        print(f"#OMIM IDs for {disease}: {len(omim_search_results)}")
        Phenotype_search_results[disease]["OMIM"]["search_params"] = Phenotype_search_parameters_OMIM[disease]
        Phenotype_search_results[disease]["OMIM"]["search_results"] = omim_search_results
        
        # print(f"Searching for {disease} in HPO...")
        # hpo_ids = search_hpo(searchTerms["hpo_st"])
        # print(f"#HPO IDs for {disease}: {len(hpo_ids)}")
        # diseases[disease]["hpo_ids"] = hpo_ids
        
        print(f"Searching for {disease} in MONDO...")
        mondo_search_results = find_disease_Phenotypes_MONDO(Phenotype_search_parameters_MONDO[disease])
        print(f"#MONDO IDs for {disease}: {len(mondo_search_results)}")
        Phenotype_search_results[disease]["MONDO"]["search_params"] = Phenotype_search_parameters_MONDO[disease]
        Phenotype_search_results[disease]["MONDO"]["search_results"] = mondo_search_results
        

        print(f"Searching for {disease} in MedGen...")
        medgen_search_results = find_disease_Phenotypes_MedGen(Phenotype_search_parameters_MedGen[disease])
        print(f"#MedGen IDs for {disease}: {len(medgen_search_results)}")
        Phenotype_search_results[disease]["MedGen"]["search_params"] = Phenotype_search_parameters_MedGen[disease]
        Phenotype_search_results[disease]["MedGen"]["search_results"] = medgen_search_results

        with open(f'{PHENOTYPES_DIR}/Disease_Phenotypes.json', 'w') as f:
            json.dump(Phenotype_search_results, f, indent=2)
    
    ClinVar_file_path = Path(ClinVar_master_file)
    if ClinVar_file_path.exists():
        ClinVar_all = pd.read_csv(ClinVar_master_file, sep="\t")
        Pheno_IDs_Ont = {ont: set() for ont in Ontologies}
        for pheno_ids_str in ClinVar_all["PhenotypeIDS"].tolist():
            _, pheno_ids_ont = extractPhenoIDs(pheno_ids_str, Ontologies=Ontologies)
            for ont in Ontologies:
                Pheno_IDs_Ont[ont].update(set(pheno_ids_ont[ont]))
        clinVar_MedGen_master = get_clinVar_MedGen_master(list(Pheno_IDs_Ont["MedGen"]))
        clinVar_OMIM_master = get_clinVar_OMIM_master(list(Pheno_IDs_Ont["OMIM"]))
        clinVar_MONDO_master = get_clinVar_MONDO_master(list(Pheno_IDs_Ont["MONDO"]))
        clinVar_master = {**clinVar_MedGen_master, **clinVar_OMIM_master, **clinVar_MONDO_master}
        with open(ClinVar_phenotypes_master_file, 'w') as f:
            json.dump(clinVar_master, f, indent=2)




# def cuid2medgen(cui):
#     # url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term={cui}&retmode=json"
#     # response = requests.get(url)
#     # response.raise_for_status()
#     # data3 = response.json()
#     # mids = data3["esearchresult"]["idlist"]
#     params = {
#         "db": "medgen",
#         "term": cui,
#         "retmode": "json",
#         "api_key": "1b42f5985b629bc2dd1cb5016ccae752db09"
#     }
#     extracted = False
#     while not extracted:
#         r = requests.get(
#             "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
#             params=params
#         )
#         if r.status_code == 429:
#             print("Rate limit hit. Sleeping...")
#             #time.sleep(2)
#             time.sleep(0.12)
#             continue
#         extracted = True
#         data = r.json()
#         mids = ["MedGen:" + mid for mid in data["esearchresult"]["idlist"]]
#         # ~8 requests/sec (safe)
#     return mids

# def search_medgenold(cuid):
#     cuids = cuid2medgen(cuid)
#     cuids = cuids+["MedGen:" + cuid]
#     API_KEY = "a76354ec-df97-4be1-82da-3e3eadd09f46"
#     BASE_URL = "https://uts-ws.nlm.nih.gov/rest/content/current/CUI"
#     url = f"{BASE_URL}/{cuid}/relations"
#     page_num = 1
#     params = {
#         "apiKey": API_KEY,
#         "pageSize": 100,
#         "pageNumber": page_num
#     }
   
#     while True:
#         params["pageNumber"] = page_num
#         response = requests.get(url, params=params)
#         response.raise_for_status()
#         data = response.json()

#         # Filter for children (narrower concepts)
#         children = [ r for r in data["result"] if r["relationLabel"] == "RB"]
    
#         for child in children:
#             #url2 = f"https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{child['ui']}/atoms?apiKey={API_KEY}"
#             rsrc = child['rootSource']
#             id = child['relatedId'].split("/")[-1]
#             url2 = f"https://uts-ws.nlm.nih.gov/rest/search/current?string={id}&sabs={rsrc}&searchType=exact&inputType=sourceUi&apiKey={API_KEY}"
#             response2 = requests.get(url2)
#             response2.raise_for_status()
#             data2 = response2.json()
#             atoms = data2["result"]["results"]
#             for atom in atoms:
#                 cuid1 = atom["ui"]
#                 cuids = cuids+["MedGen:" + cuid1]
#                     # #url3 = f"https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{cui}/atoms?apiKey={API_KEY}&pageSize=200&pageNumber={page_num}"
#                     # url3 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term={cui}&retmode=json"
#                     # response3 = requests.get(url3)
#                     # response3.raise_for_status()
#                     # data3 = response3.json()
#                     # mid = data3["esearchresult"]["idlist"]
#                 more_cuids = cuid2medgen(cuid1)
#                 cuids = cuids + more_cuids
#                 #if len(more_cuids) > 0:
#                 print(atom["ui"], atom["name"], more_cuids)
#                     # for s in sources:
#                     #     if s["rootSource"].lower() == "medgen":
#                     #         print(child["relatedId"], child["relatedIdName"])
#                     #         print(s["rootSource"], s["ui"], s["term"])
#                     #         cuids = cuids + [s["ui"]]
#         page_num += 1
#         if page_num > data["pageCount"]:
#             break    
#     cuids = list(set(cuids))
#     return cuids


#def search_mondo(mondo_id):
    # url = f"https://api.monarchinitiative.org/v3/api/ontology/descendants/{mondo_id}"
    # #url = f"https://api.monarchinitiative.org/v3/api/ontology/term/{mondo_id}/descendants"
    # response = requests.get(url)
    # response.raise_for_status()
    # data = response.json()
    # descendants = []
    # # Structure may vary slightly by API version
    # for term in data.get("terms", []):
    #     descendants.append((term["id"], term["name"]))
    # return descendants


# def search_hpo(disease_ids):
#     all_phenotypes = set()
#     #cardio_omims = {"115197", "613765", "601494"}  # add all relevant OMIM IDs

#     with open("phenotype.hpoa") as f:
#         reader = csv.reader(f, delimiter='\t')
#         for row in reader:
#             if row[0] in disease_ids:
#                 all_phenotypes.add((row[3], row[4]))  # HP ID, HP Name
#     for hp_id, name in sorted(all_phenotypes):
#         print(hp_id, "-", name)
#     IDs = [hp_id for hp_id, name in all_phenotypes]
#     return IDs



# Example: Replace with your MEDGEN CUI
# parent_cui = "C0020538"  # Example CUI
# children_terms = get_children(parent_cui)

# for child in children_terms:
#     print(child["ui"], "-", child["name"])
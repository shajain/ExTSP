from collections import defaultdict
from  ExTSP.config import DATA_DIR, DISEASE_SEMANTIC_TYPES
import gzip
from collections import defaultdict
import xml.etree.ElementTree as ET
#from Bio import Entrez
from ExTSP.config import NCBI_API_KEY
from ExTSP.config import Diseases
from ExTSP.Phenotypes.common_functions import find_in_disease_Phenotypes, find_in_ClinVar_master, get_response



raw_folder = f"{DATA_DIR}/raw"
conso_file = f"{raw_folder}/MGCONSO.RRF.gz"
rel_file = f"{raw_folder}/MGREL.RRF.gz"
#filtered_cuids_file = f"{DATA_DIR}/Phenotypes/cuids_filtered.json"
clinvar_phenotypes_file = f"{DATA_DIR}/phenotypes/clinvar_MedGen_master.json"


# ---------- Step 1: Load concept names ----------
def load_concepts():
    cui_to_name = {}
    term_to_cuis = defaultdict(set)
    #conso_file = f"{config.DATA_DIR}/MGCONSO.RRF.gz"
    with gzip.open(conso_file, 'rt', encoding='utf-8') as f:
        for i, line in enumerate(f):
            fields = line.strip().split('|')
            if i == 0:
                continue
            cui = fields[0]
            term = fields[11]  # STR field
            
            # Keep first encountered name as representative
            if cui not in cui_to_name:
                cui_to_name[cui] = term

            term_to_cuis[term.lower()].add(cui)

    return cui_to_name, term_to_cuis


# ---------- Step 2: Load relationships ----------
def load_relationships():
    parents = defaultdict(set)
    children = defaultdict(set)
    with gzip.open(rel_file, 'rt', encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split('|')
            cui1 = fields[0]
            cui2 = fields[4]
            rel = fields[3]

            # MedGen hierarchy usually encoded as:
            # PAR = parent
            # CHD = child

            if rel == "PAR":
                parents[cui1].add(cui2)
                children[cui2].add(cui1)

            elif rel == "CHD":
                children[cui1].add(cui2)
                parents[cui2].add(cui1)

    return parents, children

def get_ancestors(cuids, parents):
    ancestors = set()
    stack = list(set(cuids))
    while stack:
        current = stack.pop()
        for p in parents.get(current, []):
            if p not in ancestors:
                ancestors.add(p)
                stack.append(p)
    return ancestors

def get_descendants(cuids, children):
    descendants = set()
    stack = list(set(cuids))
    while stack:
        current = stack.pop()
        for c in children.get(current, []):
            if c not in descendants:
                descendants.add(c)
                stack.append(c)
    return descendants


def find_primary_disease_Phenotypes_MedGen(term, max_results=None):
    retmax = 100
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "medgen",
        "term": term,
        "retmode": "json",
        "retmax": retmax,
        "retstart": 0,
        "apiKey": NCBI_API_KEY
    }
    search = True
    all_uids = []
    while search:
        response = get_response(url, params=params)
        data = response.json()
        uids = data['esearchresult'].get('idlist', [])
        all_uids.extend(uids)
        #print(uids)
        if len(uids) < retmax:
            break  # No more results
        params["retstart"] += retmax
        if max_results and len(all_uids) >= max_results:
            break
    return all_uids

# ---------- Step 2: Retrieve concept details ----------
def uids2results(uids):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    batch_size = 200
    params = {
        "db": "medgen",
        "id": ",".join(uids),
        "retmode": "json",
        "apiKey": NCBI_API_KEY
    }
   
    results_dct = {}
    for i in range(0, len(uids), batch_size):
        batch_uids = uids[i:i+batch_size]
        params["id"] = ",".join(batch_uids)
        response = get_response(url, params=params)
        results = response.json()["result"]
        results = {v["conceptid"]: v for k,v in results.items() if k in batch_uids}
        results_dct.update(results)
    for cuid, entries in results_dct.items():
        root = ET.fromstring(f"<root>{entries['conceptmeta']}</root>")
        synonyms = list(set([elem.text.lower() for elem in root.findall(".//Name")]))
        results_dct[cuid]["alternative_titles"] = "|".join(synonyms)
    return results_dct
        
def filterBySemanticType(results):
    results = {cuid: entries for cuid, entries in results.items() if entries["semantictype"] and entries["semantictype"]["value"].lower() in DISEASE_SEMANTIC_TYPES}
    return results

def filterByTerm(results, filter_terms):
    summaries_neg =[]
    results1 = {}
    for cuid, attributes in results.items():
        title = attributes["title"].lower()
        description = attributes["definition"]["value"].lower() if attributes["definition"] and attributes["definition"]["value"] else ""
        if any([term.lower() in f"{title}|{description}" for term in filter_terms]):
            results1[cuid] = attributes
        else:
            # root = ET.fromstring(f"<root>{attributes['conceptmeta']}</root>")
            # synonyms = set([elem.text.lower() for elem in root.findall(".//Name")])
            # if any([term.lower() in synonym for synonym in synonyms for term in necessary_terms]):
            #     results1[cuid] = attributes
            # else:
            #     summaries_neg.append(attributes["title"])
            if any([term.lower() in attributes["alternative_titles"].lower() for term in filter_terms]):
                results1[cuid] = attributes
            else:                
                summaries_neg.append(attributes["title"])
    #print("Summaries without term:", len(summaries_neg), summaries_neg)
    return results1
        #results1 = [v for v in results if any([term.lower() in v["title"].lower() for term in necessary_terms])]
        #summaries_neg = [v["title"] for v in results if not any([term.lower() in v["title"].lower() for term in allowed_term])]

def filterResults(results, filter_terms):
    # url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    # batch_size = 200
    # params = {
    #     "db": "medgen",
    #     "id": ",".join(uids),
    #     "retmode": "json",
    #     "apiKey": NCBI_API_KEY
    # }
    # cuidss = []
    # for i in range(0, len(uids), batch_size):
    #     batch_uids = uids[i:i+batch_size]
    #     params["id"] = ",".join(batch_uids)
    #     response = get_response(url, params=params)
    #     results = response.json()["result"]
    #     results = [v for k,v in results.items() if k in batch_uids]
    #     results = [v for v in results if v["semantictype"]["value"].lower() in DISEASE_SEMANTIC_TYPES]
    results1 = filterBySemanticType(results)
        # #summaries_neg = [v["title"] for v in results if not any([term.lower() in v["title"].lower() for term in necessary_terms])]
        # summaries_neg =[]
        # results1 = []
        # for r in results:
        #     title = r["title"].lower()
        #     if any([term.lower() in title for term in necessary_terms]):
        #         results1.append(r)
        #     else:
        #         root = ET.fromstring(f"<root>{r['conceptmeta']}</root>")
        #         synonyms = set([elem.text.lower() for elem in root.findall(".//Name")])
        #         if any([term.lower() in synonym for synonym in synonyms for term in necessary_terms]):
        #             results1.append(r)
        #         else:
        #             summaries_neg.append(r["title"])
        #results1 = [v for v in results if any([term.lower() in v["title"].lower() for term in necessary_terms])]
        #summaries_neg = [v["title"] for v in results if not any([term.lower() in v["title"].lower() for term in allowed_term])]
        #print("Summaries without term:", len(summaries_neg), summaries_neg)
    results2 = filterByTerm(results1, filter_terms)
        # cuids = [v["conceptid"] for v in results2]
        # cuidss += cuids
    return results2
       
def find_disease_Phenotypes_MedGen(search_params):
    # disease_terms = {"CM": {"main": "cardiomyopathy", "necessary_terms": ["cardiomyopathy", "ventricular", "myocardial", 
    #                "arrhythmia", "heart failure", "cardiomegaly", "cardiac"]},
    #                "PCD": {"main": "primary ciliary dyskinesia", "necessary_terms": ["ciliary", "cilia"]},
    #                "PKD": {"main": "polycystic kidney disease", "necessary_terms": ["kidney", "renal", "polycystic" ,"nephro"]},
    #                "PH": {"main": "pulmonary hypertension", "necessary_terms": ["pulmonary", "hypertension", "lung", 
    #                                                                             "arterial", "pulmonic", "right heart"]}
    #                                                                             }

    filter_terms = search_params["filter_terms"]
    main_term = search_params["search_str"]
    uids = find_primary_disease_Phenotypes_MedGen(main_term, max_results=None)
    results = uids2results(uids)
    results_filtered = filterResults(results, filter_terms)
    results_filtered = {c:{**entries, "hierarchy": "main"} for c, entries in results_filtered.items()}
    cuids_entry = [c for c, entries in results_filtered.items()]
    #cuids = uids2cuids(uids, necessary_terms=necessary_terms)
    cuids = list(results_filtered.keys())
    parents_relationships, children_relationships = load_relationships()
    descendants = get_descendants(cuids, children_relationships)
    descendants = list(set(descendants)-set(cuids))
    ancestors = get_ancestors(cuids, parents_relationships)
    ancestors = list(set(ancestors)-set(cuids))
    ancestors = list(set(ancestors)-set(descendants))
    descendants_filtered = filterBySemanticType(uids2results(cuids2uids(list(descendants))))
    ancestors_filtered = filterBySemanticType(uids2results(cuids2uids(list(ancestors))))
    descendants_filtered = {c:{**entries, "hierarchy": "descendant"} for c, entries in descendants_filtered.items()}
    ancestors_filtered = {c:{**entries, "hierarchy": "ancestor"} for c, entries in ancestors_filtered.items()}
    all_results = {**results_filtered, **descendants_filtered, **ancestors_filtered}
    all_results = {f"MedGen:{c}":{"uid": entries["uid"], "title": entries["title"], 
                                  "semantictype": entries["semantictype"]["value"], 
                                  "hierarchy": entries["hierarchy"],
                                  "alternative_titles": entries["alternative_titles"]} for c, entries in all_results.items()}
    # cuids_all = set(cuids + list(descendants))
    # cuids_all = ["MedGen:" + cuid for cuid in cuids_all]
    return all_results


def extractRelevantFields(result_value):
    important_fields = ["uid", "title", "semantictype", "hierarchy", "alternative_titles"] 
    result_value1 = {field: result_value.get(field, None) for field in important_fields} 
    for field in important_fields:
        if field == "semantictype" and "value" in result_value1[field]:
            result_value1[field] = result_value1[field]["value"]
        # elif field in result_value:
        #     result_value1[field] = result_value[field]
    return result_value1


    

def MedGen_extraInfo(cuid):
    info = find_in_ClinVar_master(cuid, "MedGen")
    if info:
        return info
    for disease in Diseases.keys():
        info = find_in_disease_Phenotypes(cuid, disease)
        if info:
            return info
    uids = cuids2uids([cuid])
    if uids:
        results = uids2results(uids)
        if results:
            results = {c: extractRelevantFields(entries) for c, entries in results.items()}
            return results[cuid]
    return None


def cuids2uids(cuids):
    #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term={query}&retmode=json"
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # Build OR query
    
    batch_size = 200
    params = {
        "db": "medgen",
        # "term": query,
        "retmode": "json",
        "apiKey": NCBI_API_KEY
    }
    uidss = []
    for cuid in cuids:
        if cuid.lower().startswith("medgen:"):
            cuid = cuid.split(":")[1]
        params["term"] = cuid
        response = get_response(url, params=params)
        uids = response.json()["esearchresult"]["idlist"]
        if uids:
            uid = uids[0]
        else:
            uid = None
        uidss.append(uid)
    return uidss


def get_clinVar_MedGen_master(cuids):
    uids = cuids2uids(cuids)
    cuids1 = [cuid for cuid, uid in zip(cuids, uids) if uid is not None]
    uids1 = [uid for uid in uids if uid is not None]
    print(f"Filtered out {len(cuids) - len(cuids1)} cuids without corresponding uids")
    uids1 = list(set(uids1))
    results = uids2results(uids1)
    assert len(results) == len(uids1), "Mismatch between results and uids"
    results_filtered = filterBySemanticType(results)
    print(f"Filtered out {len(results) - len(results_filtered)} results that do not correspond to diseases")
    cuids_filtered = {f"MedGen:{id}": extractRelevantFields(fields) for id, fields in results_filtered.items()}
    # with open(clinvar_phenotypes_file, "w") as f:
    #         json.dump(cuids_filtered, f)
    #         print(f"Saved {len(cuids_filtered)} filtered cuids to file")
    return cuids_filtered

# def get_clinVar_MedGen_master(cuids):
#     file_path = Path(clinvar_phenotypes_file)
#     if file_path.exists():
#         with open(clinvar_phenotypes_file, "r") as f:
#             cuids_filtered = json.load(f)
#             print(f"Loaded {len(cuids_filtered)} filtered cuids from file")
#     else:
#         cuids_filtered = clinVar_MedGen_master(cuids)
#     return cuids_filtered

# def filter_clinVar_MedGenIDs(cuids):
#     file_path = Path(clinvar_phenotypes_file)
#     if file_path.exists():
#         with open(clinvar_phenotypes_file, "r") as f:
#             cuids_filtered = json.load(f)
#             print(f"Loaded {len(cuids_filtered)} filtered cuids from file")
#     else:
#         cuids_filtered = create_clinVar_MedGen_master(cuids)
#     return list(cuids_filtered.keys())
    

#---------- Run ----------
if __name__ == "__main__":
    # necessary_terms = {"cardiomyopathy":["cardiomyopathy", "ventricular", "myocardial", 
    #                "arrhythmia", "heart failure", "cardiomegaly", "cardiac"]}
    all_disease_phenotypes = find_disease_Phenotypes_MedGen("CM")
    print(all_disease_phenotypes)
    




# def cuids2uids(cuids):
#     #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term={query}&retmode=json"
#     url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#     # Build OR query
    
#     batch_size = 200
#     params = {
#         "db": "medgen",
#         # "term": query,
#         "retmode": "json",
#         "apiKey": NCBI_API_KEY
#     }
#     uidss = []
#     for i in range(0, len(cuids), batch_size):
#         batch_cuids = cuids[i:i+batch_size]
#         query = " OR ".join([f"{cid}" for cid in batch_cuids])
#         params["term"] = query
#         r = requests.get(url, params=params)
#         try:
#             r.raise_for_status()
#         except requests.exceptions.HTTPError as e:
#             if wait_for_rate_limit(r):
#                 continue
#         results = r.json()["result"]
#         uids = results["idlist"]
#         uidss += uids
#     return uidss


#def find_in_Stored_Phenotypes(cuid, file_path):
#     if "medgen" not in cuid.lower():
#         cuid = "MedGen:" + cuid
#     if file_path.exists():
#         with open(file_path, "r") as f:
#             phenotypes = json.load(f)
#             if cuid in phenotypes:
#                 return phenotypes[cuid]
#     return None

# def find_in_ClinVar_master(cuid):
#     file_path = Path(clinvar_phenotypes_file)
#     return find_in_Stored_Phenotypes(cuid, file_path)


# def find_in_disease_Phenotypes(cuid, disease):
#     disease_phenotype_file = f"{DATA_DIR}/Phenotypes/{disease}_Phenotypes.json"
#     file_path = Path(disease_phenotype_file)
#     return find_in_Stored_Phenotypes(cuid, file_path)
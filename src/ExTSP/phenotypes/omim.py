from ExTSP.config import DATA_DIR, OMIM_API_KEY
from ExTSP.Phenotypes.common_functions import get_response
Phenotypes_folder = f"{DATA_DIR}/phenotypes"
clinvar_phenotypes_file = f"{Phenotypes_folder}/clinvar_OMIM_master.json"

def find_disease_Phenotypes_OMIM(search_params, start=0, limit=100):
    BASE_URL = "https://api.omim.org/api"
    params = {
        "search": f"title:{search_params['search_str']}",
        "include": "title",
        "format": "json",
        "apiKey": OMIM_API_KEY,
        "start": start,
        "limit": limit,  
    }
    params2 = {k:v for k,v in params.items() if k!="search"}
    IDs = {}
    while True:
        params["start"] = start
        response = get_response(f"{BASE_URL}/entry/search", params=params)
        if response.status_code == 200:
            data = response.json()
        else:
            break
        entries = data['omim']['searchResponse']['entryList']
        params2["mimNumber"] = ",".join([e['entry']['mimNumber'] for e in entries])
        response2 = get_response(f"{BASE_URL}/entry", params=params2)
        if response2.status_code == 200:
            data2 = response2.json()
        else:
            break
        entries2 = data2['omim']['entryList']
        assert len(entries) == len(entries2), "Mismatch in number of entries returned by search and entry endpoints"
        
        for entry, entry2 in zip(entries, entries2):
            e = entry['entry']
            e2 = entry2['entry']
            #print(e['mimNumber'], e['titles']['preferredTitle'])
            if "alternativeTitles" in e2["titles"]:
                #print("ho")  # filter out very long titles that are unlikely to be relevant
                IDs['OMIM:'+str(e['mimNumber'])] = {"title": e['titles']['preferredTitle'], "alternative_titles": e2['titles']['alternativeTitles']}
            else:
                IDs['OMIM:'+str(e['mimNumber'])] = {"title": e['titles']['preferredTitle'], "alternative_titles": ""}
        start += len(entries)
        if start >= data['omim']['searchResponse']['totalResults']:
            break
    return IDs




def find_OMIM_id_api(omim_ids):
    omim_ids = [omim_ids] if isinstance(omim_ids, str) else omim_ids
    omim_ids = [id.split(":")[1] if id.startswith("OMIM:") else id for id in omim_ids]
    BASE_URL = "https://api.omim.org/api"
    batch_size = 100
    params = {
        "include": "title",
        "format": "json",
        "apiKey": OMIM_API_KEY,
    }
    omim_ids_info = {}
    for i in range(0, len(omim_ids), batch_size):
        batch_ids = omim_ids[i:i+batch_size]
        params["mimNumber"] = ",".join(batch_ids)
        response = get_response(f"{BASE_URL}/entry", params=params)
        if response.status_code == 200:
            data = response.json()
            entries = data['omim']['entryList']
        else:
            raise Exception(f"Error: {response.status_code}")
            return None
        for entry in entries:
            omim_ids_info[entry['entry']['mimNumber']] = {"title": entry['entry']['titles']['preferredTitle'], "alternative_titles": entry['entry']['titles'].get('alternativeTitles', "")}
    return omim_ids_info

def OMIM_extraInfo(omim_id, use_api_only=False):
    if not use_api_only:
        info = find_in_ClinVar_master(omim_id, "OMIM")
        if info:
            return info
        info = find_in_disease_Phenotypes(omim_id, ontology="OMIM")
        if info:
            return info
    info = find_OMIM_id_api(omim_id)
    return info
   
    
def get_clinVar_OMIM_master(omim_ids):
    omim_ids_info = {id: OMIM_extraInfo(id, use_api_only=True) for id in omim_ids}
    omim_ids_info_filtered = {id: info for id, info in omim_ids_info.items() if info}
    # with open(clinvar_phenotypes_file, "w") as f:
    #         json.dump(omim_ids_info_filtered, f)
    #         print(f"Saved {len(omim_ids_info_filtered)} filtered OMIM IDs to file")
    return omim_ids_info_filtered



if __name__ == "__main__":
    Ids = find_disease_Phenotypes_OMIM("cardiomyopathy")
    print(Ids)




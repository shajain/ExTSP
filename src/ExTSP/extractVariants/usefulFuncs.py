
import re
from ExTSP.config import Ontologies

def extractPhenoIDs(pheno_ids_str, Ontologies=Ontologies):
    #pid_str = pid_str.replace("Human Phenotype Ontology", "HP")
    pheno_ids_str = "".join(pheno_ids_str.split())
    #pid_str = pid_str.replace("HP:HP", "HP")
    for ont in Ontologies:
        pheno_ids_str = pheno_ids_str.replace(f"{ont}:{ont}", f"{ont}")
    # pid_str = pid_str.replace("OMIM:OMIM", "OMIM")
    # pid_str = pid_str.replace("MONDO:MONDO", "MONDO")
    # pid_str = pid_str.replace("MedGen:MedGen", "MedGen")
    pheno_ids = re.split(r"[\,\;\|]", pheno_ids_str) 
    # pids = [pid.strip() for pid in pids if pid.strip()]
    pheno_ids_ont = {}
    for ont in Ontologies:
        pheno_ids_ont[ont] = [pid for pid in pheno_ids if pid.lower().startswith(ont.lower())]
    pheno_ids = [pheno_id for ont in Ontologies for pheno_id in pheno_ids_ont[ont]]
    return pheno_ids, pheno_ids_ont



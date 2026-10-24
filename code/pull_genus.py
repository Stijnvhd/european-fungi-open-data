import requests, json, time, os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
K=json.load(open("data/keys.json"))
G0={"REF":[186],
   "SENS":[d["key"] for d in K["GEN_SENS"]],
   "TOL":[d["key"] for d in K["GEN_TOL"]]+[d["key"] for d in K["SP_TOL"]],
   "CANTH":[5249504],
   "ECM_ALL":[d["key"] for d in K["GEN_SENS"]]+[d["key"] for d in K["GEN_TOL"]]+[d["key"] for d in K["SP_TOL"]]}
G={("G_"+d["name"]):[d["key"]] for g in ["GEN_SENS","GEN_TOL","SP_TOL"] for d in K[g]}
C=["NL","BE","DK","GB","DE","NO","SE","FI","FR"]
S=requests.Session()
import sys
SRC=sys.argv[1]
CORE = {"DK": ["84d26682-f762-11e1-a439-00145eb45e9a"], "GB": ["b3604297-f438-4b5c-b654-7976a909f24e", "823cb2e1-562a-433c-be9d-d18e44d09a64"],
        "BE": ["039f9753-1e6a-4ffc-83f9-daac3686d165"], "SE": ["38b4c89f-584c-41bb-bd8f-cd1def33e92f"], "NO": ["b124e1e0-4755-430f-9eab-894f25a9b59c"],
        "NL": ["8a863029-f435-446a-821e-275f4f641165"], "DE": ["6ac3f774-d9fb-4796-b3e9-92bf6c81c084"], "FI": ["2bc4c2db-1dfa-419c-93ce-6602c5ae8c99"],
        "FR": ["bd8ee748-4ae8-4956-98d9-c1b8e8b2b9c9"]}
def q(grp,c,y):
    p=[("taxonKey",k) for k in G[grp]]+[("country",c),("year",y),("limit",0),("facet","gadmLevel2Gid"),("facetLimit",2000),("hasCoordinate","true"),("hasGeospatialIssue","false"),("occurrenceStatus","PRESENT")]+[("basisOfRecord",b_) for b_ in ["HUMAN_OBSERVATION","PRESERVED_SPECIMEN","OCCURRENCE","OBSERVATION"]]+([("datasetKey",k_) for k_ in CORE[c]] if SRC=="core" else [])
    for a in range(8):
        try:
            r=S.get("https://api.gbif.org/v1/occurrence/search",params=p,timeout=90)
            if r.status_code==200:
                j=r.json(); return [(grp,c,y,x["name"],x["count"]) for x in (j["facets"][0]["counts"] if j["facets"] else [])]
        except Exception: pass
        time.sleep(3*(a+1))
    print("FAIL",grp,c,y,flush=True); return []
jobs=[(g,c,y) for g in G for c in C for y in range(1985,2026)]
rows=[]
with ThreadPoolExecutor(4) as ex:
    for i,res in enumerate(ex.map(lambda a:q(*a),jobs)):
        rows+=res
        if i%100==0: print(i,len(jobs),flush=True)
pd.DataFrame(rows,columns=["grp","country","year","gid","n"]).to_parquet(f"data/ngroups_genus_{SRC}.parquet");print("done")

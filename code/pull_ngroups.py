import requests, json, time, os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
K=json.load(open("data/keys.json"))
G={"REF":[186],
   "SENS":[d["key"] for d in K["GEN_SENS"]],
   "TOL":[d["key"] for d in K["GEN_TOL"]]+[d["key"] for d in K["SP_TOL"]],
   "CANTH":[5249504],
   "ECM_ALL":[d["key"] for d in K["GEN_SENS"]]+[d["key"] for d in K["GEN_TOL"]]+[d["key"] for d in K["SP_TOL"]]}
C=["NL","BE","DK","GB","DE","NO","SE","FI","FR"]
S=requests.Session()
def q(grp,c,y):
    p=[("taxonKey",k) for k in G[grp]]+[("country",c),("year",y),("limit",0),("facet","gadmLevel2Gid"),("facetLimit",2000),("hasCoordinate","true"),("hasGeospatialIssue","false"),("occurrenceStatus","PRESENT")]
    for a in range(8):
        try:
            r=S.get("https://api.gbif.org/v1/occurrence/search",params=p,timeout=90)
            if r.status_code==200:
                j=r.json(); return [(grp,c,y,x["name"],x["count"]) for x in (j["facets"][0]["counts"] if j["facets"] else [])]
        except Exception: pass
        time.sleep(3*(a+1))
    print("FAIL",grp,c,y,flush=True); return []
jobs=[(g,c,y) for g in G for c in C for y in range(1980,2026)]
rows=[]
with ThreadPoolExecutor(4) as ex:
    for i,res in enumerate(ex.map(lambda a:q(*a),jobs)):
        rows+=res
        if i%100==0: print(i,len(jobs),flush=True)
pd.DataFrame(rows,columns=["grp","country","year","gid","n"]).to_parquet("data/ngroups_gadm2.parquet");print("done")

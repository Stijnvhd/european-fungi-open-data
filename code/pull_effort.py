import requests, time
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
REFS={"PLANTS":7707728,"AGARICO":186,"BIRDS":212}
C=["NL","BE","DK","GB","DE","NO","SE","FI","FR"]
S=requests.Session()
def q(ref,c,y,m):
    p={"taxonKey":REFS[ref],"country":c,"year":y,"month":m,"limit":0,"facet":"day","facetLimit":40,"hasCoordinate":"true","occurrenceStatus":"PRESENT"}
    for a in range(8):
        try:
            r=S.get("https://api.gbif.org/v1/occurrence/search",params=p,timeout=90)
            if r.status_code==200:
                j=r.json(); return [(ref,c,y,m,int(x["name"]),x["count"]) for x in (j["facets"][0]["counts"] if j["facets"] else [])]
        except Exception: pass
        time.sleep(3*(a+1))
    print("FAIL",ref,c,y,m,flush=True); return []
jobs=[(r,c,y,m) for r in ["PLANTS","AGARICO"] for c in C for y in range(1980,2026) for m in range(1,13)]
rows=[]
with ThreadPoolExecutor(4) as ex:
    for i,res in enumerate(ex.map(lambda a:q(*a),jobs)):
        rows+=res
        if i%500==0: print(i,len(jobs),flush=True)
pd.DataFrame(rows,columns=["ref","country","year","month","day","n"]).to_parquet("data/effort_daily.parquet");print("done")

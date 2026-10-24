import requests, json, os, time, sys
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
K=json.load(open("data/keys.json"))
C=["NL","BE","DK","GB","NO","SE","DE","FI","FR"]
OUT="data/records"; os.makedirs(OUT,exist_ok=True)
S=requests.Session()
def get(params):
    for a in range(8):
        try:
            r=S.get("https://api.gbif.org/v1/occurrence/search",params=params,timeout=60)
            if r.status_code==200: return r.json()
        except Exception as e: pass
        time.sleep(2*(a+1))
    raise RuntimeError(params)
FIELDS=["key","year","month","day","eventDate","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","datasetKey","basisOfRecord","countryCode","recordedBy","institutionCode"]
def pull_range(key,c,y0,y1):
    base={"taxonKey":key,"country":c,"year":f"{y0},{y1}","hasCoordinate":"true","hasGeospatialIssue":"false","occurrenceStatus":"PRESENT","limit":300}
    n=get({**base,"limit":0})["count"]
    if n==0: return []
    if n>99000 and y1>y0:
        m=(y0+y1)//2
        return pull_range(key,c,y0,m)+pull_range(key,c,m+1,y1)
    rows=[];off=0
    while off<n and off<100000:
        j=get({**base,"offset":off})
        for x in j["results"]:
            rows.append({f:x.get(f) for f in FIELDS})
        if j.get("endOfRecords"): break
        off+=300
    return rows
def job(g,d,c):
    fn=f"{OUT}/{d['key']}_{c}.parquet"
    if os.path.exists(fn): return fn,0
    rows=pull_range(d["key"],c,1970,2025)
    df=pd.DataFrame(rows,columns=FIELDS); df["taxon"]=d["name"]; df["guild"]=g
    for col in ["recordedBy","institutionCode","eventDate","basisOfRecord","datasetKey","countryCode"]: df[col]=df[col].astype("string")
    df.to_parquet(fn); return fn,len(df)
jobs=[(g,d,c) for g in ["ECM","SAP"] for d in K[g] for c in C]
t0=time.time();done=0
with ThreadPoolExecutor(8) as ex:
    fs=[ex.submit(job,*j) for j in jobs]
    for f in as_completed(fs):
        done+=1
        try: fn,n=f.result(); print(done,len(jobs),fn,n,round(time.time()-t0),flush=True)
        except Exception as e: print("ERR",e,flush=True)

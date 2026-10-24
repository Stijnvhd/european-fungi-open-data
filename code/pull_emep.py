from netCDF4 import Dataset
import numpy as np, pandas as pd, time, os
B="https://thredds.met.no/thredds/dodsC/data/EMEP/2025_Reporting/"
rows=[]
for y in range(1990,2024):
    fn=f"EMEP01_rv5.6_year.{y}met_{y}emis_rep2025.nc" if y<2023 else "EMEP01_rv5.6_year.2023met_2023emis.nc"
    for a in range(5):
        try:
            d=Dataset(B+fn); break
        except Exception as e: print("retry",y,e); time.sleep(10)
    la=d["lat"][:];lo=d["lon"][:]
    i0,i1=np.searchsorted(la,44.0),np.searchsorted(la,72.0)
    j0,j1=np.searchsorted(lo,-11.0),np.searchsorted(lo,32.0)
    tot=None;parts={}
    for v in ["DDEP_OXN_m2Grid","WDEP_OXN","DDEP_RDN_m2Grid","WDEP_RDN"]:
        a=np.asarray(d[v][0,i0:i1,j0:j1],dtype=float); parts[v]=a
    LA,LO=np.meshgrid(la[i0:i1],lo[j0:j1],indexing="ij")
    df=pd.DataFrame({"lat":LA.ravel(),"lon":LO.ravel(),**{k:v.ravel() for k,v in parts.items()}})
    df["cell_lat"]=np.floor(df.lat/0.5)*0.5;df["cell_lon"]=np.floor(df.lon/0.5)*0.5
    g=df.groupby(["cell_lat","cell_lon"])[list(parts)].mean().reset_index();g["year"]=y
    rows.append(g);print(y,len(g),flush=True);d.close()
out=pd.concat(rows)
out["oxN"]=out.DDEP_OXN_m2Grid+out.WDEP_OXN; out["redN"]=out.DDEP_RDN_m2Grid+out.WDEP_RDN
out["totN_kg_ha"]=(out.oxN+out.redN)/100.0  # mgN/m2 -> kgN/ha
out.to_parquet("data/emep_ndep_halfdeg.parquet");print("done")

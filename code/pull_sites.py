import requests, time, json, os
import pandas as pd
import numpy as np

P = pd.read_csv("data/sites.csv")
S = requests.Session()


def get(url, params, tries=6):
    for k in range(tries):
        try:
            r = S.get(url, params=params, timeout=120, headers={"Accept": "application/json"})
            if r.status_code == 200:
                return r.json()
            if r.status_code == 429:
                time.sleep(30)
        except Exception:
            pass
        time.sleep(5 * (k + 1))
    return None


# Daily climate 1981-2025 per site (NASA POWER, MERRA-2 based)
clim = []
for _, s in P.iterrows():
    j = get("https://power.larc.nasa.gov/api/temporal/daily/point",
            {"parameters": "T2M,PRECTOTCORR", "community": "AG", "longitude": round(s.lon, 3), "latitude": round(s.lat, 3),
             "start": "19810101", "end": "20251231", "format": "JSON"})
    if j is None:
        print("climate FAIL", s.gid, flush=True)
        continue
    pr = j["properties"]["parameter"]
    d = pd.DataFrame({"time": list(pr["T2M"].keys()), "t2m": list(pr["T2M"].values()), "prec": [pr["PRECTOTCORR"][k] for k in pr["T2M"]]})
    d["gid"] = s.gid
    d["country"] = s.country
    clim.append(d)
    print("climate", s.gid, len(d), flush=True)
pd.concat(clim).to_parquet("data/site_climate.parquet")

# MODIS MCD12Q2 land surface phenology 2001-2023, 5x5 km window
bands = ["MidGreenup.Num_Modes_01", "MidGreendown.Num_Modes_01", "Dormancy.Num_Modes_01"]
rows = []
for _, s in P.iterrows():
    dts = get("https://modis.ornl.gov/rst/api/v1/MCD12Q2/dates", {"latitude": s.lat, "longitude": s.lon})
    if dts is None:
        print("modis dates FAIL", s.gid, flush=True)
        continue
    md = [d["modis_date"] for d in dts["dates"]]
    for b in bands:
        for i in range(0, len(md), 10):
            ch = md[i:i + 10]
            j = get("https://modis.ornl.gov/rst/api/v1/MCD12Q2/subset",
                    {"latitude": s.lat, "longitude": s.lon, "band": b, "startDate": ch[0], "endDate": ch[-1],
                     "kmAboveBelow": 2, "kmLeftRight": 2})
            if j is None:
                print("modis FAIL", s.gid, b, ch[0], flush=True)
                continue
            for sub in j.get("subset", []):
                v = np.array(sub["data"], dtype=float)
                v = v[(v > 0) & (v < 32767)]
                if len(v) == 0:
                    continue
                yr = int(sub["calendar_date"][:4])
                doy = v - (pd.Timestamp(f"{yr}-01-01") - pd.Timestamp("1970-01-01")).days + 1
                doy = doy[(doy > 0) & (doy < 400)]
                rows.append(dict(gid=s.gid, country=s.country, band=b.split(".")[0], year=yr,
                                 doy_median=float(np.median(doy)) if len(doy) else np.nan, npix=int(len(doy))))
    print("modis", s.gid, flush=True)
pd.DataFrame(rows).to_parquet("data/site_modis.parquet")
print("done")

import requests, time
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

P = pd.read_csv("data/sites.csv")
BANDS = ["MidGreenup.Num_Modes_01", "MidGreendown.Num_Modes_01", "Dormancy.Num_Modes_01"]
BASE = "https://modis.ornl.gov/rst/api/v1/MCD12Q2/"


def get(url, params):
    for k in range(6):
        try:
            r = requests.get(url, params=params, timeout=180, headers={"Accept": "application/json"})
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(5 * (k + 1))
    return None


def site(s):
    rows = []
    dts = get(BASE + "dates", {"latitude": s.lat, "longitude": s.lon})
    if dts is None:
        print("FAIL dates", s.gid, flush=True)
        return rows
    md = [d["modis_date"] for d in dts["dates"]]
    for b in BANDS:
        for i in range(0, len(md), 10):
            ch = md[i:i + 10]
            j = get(BASE + "subset", {"latitude": s.lat, "longitude": s.lon, "band": b, "startDate": ch[0],
                                      "endDate": ch[-1], "kmAboveBelow": 2, "kmLeftRight": 2})
            if j is None:
                print("FAIL", s.gid, b, ch[0], flush=True)
                continue
            for sub in j.get("subset", []):
                v = np.array(sub["data"], dtype=float)
                v = v[(v > 0) & (v < 32767)]
                yr = int(sub["calendar_date"][:4])
                doy = v - (pd.Timestamp(f"{yr}-01-01") - pd.Timestamp("1970-01-01")).days + 1
                doy = doy[(doy > 0) & (doy < 400)]
                if len(doy):
                    rows.append(dict(gid=s.gid, country=s.country, band=b.split(".")[0], year=yr,
                                     doy_median=float(np.median(doy)), npix=int(len(doy))))
    print("done", s.gid, flush=True)
    return rows


with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(site, [s for _, s in P.iterrows()]))
pd.DataFrame([r for rr in res for r in rr]).to_parquet("data/site_modis.parquet")
print("finished")

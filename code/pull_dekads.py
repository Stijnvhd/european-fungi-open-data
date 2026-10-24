import requests, json, time, sys
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

K = json.load(open("data/keys.json"))
src = sys.argv[1]  # "fb" (all fruitbody records) or "core" (main national recording scheme only)
CORE = {"DK": ["84d26682-f762-11e1-a439-00145eb45e9a"], "GB": ["b3604297-f438-4b5c-b654-7976a909f24e", "823cb2e1-562a-433c-be9d-d18e44d09a64"],
        "BE": ["039f9753-1e6a-4ffc-83f9-daac3686d165"], "SE": ["38b4c89f-584c-41bb-bd8f-cd1def33e92f"], "NO": ["b124e1e0-4755-430f-9eab-894f25a9b59c"],
        "NL": ["8a863029-f435-446a-821e-275f4f641165"], "DE": ["6ac3f774-d9fb-4796-b3e9-92bf6c81c084"], "FI": ["2bc4c2db-1dfa-419c-93ce-6602c5ae8c99"],
        "FR": ["bd8ee748-4ae8-4956-98d9-c1b8e8b2b9c9"]}
BASIS = ["HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OCCURRENCE", "OBSERVATION"]
T = [(d["name"], g, d["key"]) for g in ["ECM", "SAP"] for d in K[g]]
T += [("Agaricomycetes", "REF", 186), ("Tracheophyta", "REF", 7707728)]
C = ["NL", "BE", "DK", "GB", "DE", "NO", "SE", "FI", "FR"]
DK = [(m, a, b, i) for m in range(1, 13) for i, (a, b) in enumerate([(1, 10), (11, 20), (21, 31)])]
S = requests.Session()


def q(t, c, dk):
    name, g, key = t
    m, a, b, i = dk
    p = [("taxonKey", key), ("country", c), ("month", m), ("day", f"{a},{b}"), ("limit", 0),
         ("facet", "year"), ("facetLimit", 150), ("hasCoordinate", "true"),
         ("hasGeospatialIssue", "false"), ("occurrenceStatus", "PRESENT")] + [("basisOfRecord", b_) for b_ in BASIS]
    if src == "core":
        p += [("datasetKey", k_) for k_ in CORE[c]]
    for k in range(8):
        try:
            r = S.get("https://api.gbif.org/v1/occurrence/search", params=p, timeout=90)
            if r.status_code == 200:
                j = r.json()
                fc = j["facets"][0]["counts"] if j["facets"] else []
                return [(name, g, c, int(x["name"]), m, i, x["count"]) for x in fc]
        except Exception:
            pass
        time.sleep(3 * (k + 1))
    print("FAIL", name, c, dk, flush=True)
    return []


jobs = [(t, c, dk) for t in T for c in C for dk in DK]
rows = []
with ThreadPoolExecutor(6) as ex:
    for i, res in enumerate(ex.map(lambda a: q(*a), jobs)):
        rows += res
        if i % 1000 == 0:
            print(i, len(jobs), flush=True)
pd.DataFrame(rows, columns=["taxon", "guild", "country", "year", "month", "dekad", "n"]).to_parquet(f"data/dekads_{src}.parquet")
print("done")

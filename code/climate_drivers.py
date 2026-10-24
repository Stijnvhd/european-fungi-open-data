"""Does the weather explain the change in the fungal season?"""
import json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

CL = pd.read_parquet("data/site_climate.parquet")
CL["date"] = pd.to_datetime(CL.time, format="%Y%m%d")
CL = CL[(CL.t2m > -900) & (CL.prec > -900)].copy()
CL["year"] = CL.date.dt.year
CL["month"] = CL.date.dt.month

def seas(months, var, how):
    x = CL[CL.month.isin(months)].groupby(["country", "gid", "year"])[var]
    x = x.mean() if how == "mean" else x.sum()
    return x.groupby(["country", "year"]).mean()

K = pd.DataFrame({
    "t_jja": seas([6, 7, 8], "t2m", "mean"),
    "t_so": seas([9, 10], "t2m", "mean"),
    "t_on": seas([10, 11], "t2m", "mean"),
    "p_ja": seas([7, 8], "prec", "sum"),
    "p_as": seas([8, 9], "prec", "sum"),
    "p_so": seas([9, 10], "prec", "sum"),
}).reset_index()
K.to_csv("out/country_climate.csv", index=False)
res = {}
# warming per decade of each season, pooled across countries
for v in ["t_jja", "t_so", "t_on", "p_ja", "p_as", "p_so"]:
    m = smf.ols(f"{v} ~ I((year-2000)/10) + C(country)", data=K[K.year.between(1981, 2025)]).fit(cov_type="HC1")
    res[f"trend_{v}"] = dict(per_decade=float(m.params.iloc[-1]), p=float(m.pvalues.iloc[-1]))

M = pd.read_parquet("out/phenology_metrics.parquet")
out = {}
for ver in ["effort", "core"]:
    d = M[(M.version == ver) & (M.year >= 1981)].merge(K, on=["country", "year"])
    d["sc"] = d.taxon + "_" + d.country
    d["w"] = np.sqrt(d.n.clip(upper=2000))
    d["yr"] = (d.year - 2000) / 10
    for gu in ["ECM", "SAP"]:
        dd = d[d.guild == gu]
        for var in ["start", "middle", "end"]:
            base = smf.wls(f"{var} ~ yr + C(sc)", data=dd, weights=dd.w).fit(cov_type="cluster", cov_kwds={"groups": dd.sc})
            full = smf.wls(f"{var} ~ yr + t_jja + t_so + t_on + p_ja + p_as + p_so + C(sc)", data=dd, weights=dd.w).fit(
                cov_type="cluster", cov_kwds={"groups": dd.sc})
            clim = smf.wls(f"{var} ~ t_jja + t_so + t_on + p_ja + p_as + p_so + C(sc)", data=dd, weights=dd.w).fit(
                cov_type="cluster", cov_kwds={"groups": dd.sc})
            out[f"{ver}|{gu}|{var}"] = dict(
                trend_alone=float(base.params["yr"]), trend_with_climate=float(full.params["yr"]),
                trend_with_climate_p=float(full.pvalues["yr"]),
                **{f"{c}": float(clim.params[c]) for c in ["t_jja", "t_so", "t_on", "p_ja", "p_as", "p_so"]},
                **{f"{c}_p": float(clim.pvalues[c]) for c in ["t_jja", "t_so", "t_on", "p_ja", "p_as", "p_so"]},
                n=int(len(dd)))
res["models"] = out
json.dump(res, open("out/climate_results.json", "w"), indent=1)
for k, v in res.items():
    if k != "models":
        print(k, v)
for k, v in out.items():
    print(k, {a: round(b, 2) for a, b in v.items()})

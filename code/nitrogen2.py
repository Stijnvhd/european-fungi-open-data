"""Nitrogen deposition and the ectomycorrhizal fungi that avoid it.

Unit of analysis: GADM level-2 administrative unit x year, 1990-2025.
Response: the share of nitrogen-sensitive genera among all ectomycorrhizal
records of the indicator genera (sensitive + tolerant). Because both halves
of the ratio are recorded by the same people on the same forays, observer
effort largely cancels. iNaturalist records are removed as a robustness
check because that platform arrived late and recruits a different recorder.
"""
import json
import numpy as np
import pandas as pd
import pyfixest as pf

import sys
SRC = sys.argv[1]
d = pd.read_parquet(f"data/ngroups_gadm2_{SRC}.parquet")
di = pd.read_parquet("data/ngroups_gadm2_inat.parquet")
cen = pd.read_parquet("data/gadm2_centroids.parquet")
dep = pd.read_parquet("data/emep_ndep_halfdeg.parquet")

def wide(x):
    w = x.pivot_table(index=["country", "gid", "year"], columns="grp", values="n", aggfunc="sum", fill_value=0)
    return w.reset_index()

W = wide(d)
Wi = wide(di).rename(columns={g: g + "_inat" for g in ["REF", "SENS", "TOL", "CANTH", "ECM_ALL"]})
W = W.merge(Wi, on=["country", "gid", "year"], how="left").fillna(0)
for g in ["REF", "SENS", "TOL", "CANTH", "ECM_ALL"]:
    W[g + "_noinat"] = (W[g] - W[g + "_inat"]).clip(lower=0)
W = W.merge(cen[["gid", "lat", "lon", "name", "area_km2", "iso"]], on="gid", how="inner")
ISO = {"NL": "NLD", "BE": "BEL", "DK": "DNK", "GB": "GBR", "DE": "DEU", "NO": "NOR", "SE": "SWE", "FI": "FIN", "FR": "FRA"}
W = W[W.iso == W.country.map(ISO)]
W["cell_lat"] = np.floor(W.lat / 0.5) * 0.5
W["cell_lon"] = np.floor(W.lon / 0.5) * 0.5

# deposition history per cell; years after 2023 carry the 2023 value
dep = dep[["cell_lat", "cell_lon", "year", "totN_kg_ha", "redN", "oxN"]].copy()
last = dep[dep.year == 2023]
ext = pd.concat([last.assign(year=y) for y in [2024, 2025]])
dep = pd.concat([dep, ext]).sort_values(["cell_lat", "cell_lon", "year"])
g = dep.groupby(["cell_lat", "cell_lon"])
dep["dep0"] = dep.totN_kg_ha
for k in [3, 5, 10]:
    dep[f"dep_ma{k}"] = g.totN_kg_ha.transform(lambda s: s.rolling(k, min_periods=1).mean())
dep["dep_lag5"] = g.totN_kg_ha.shift(5)
dep["dep_lag10"] = g.totN_kg_ha.shift(10)
dep["dep_cum"] = g.totN_kg_ha.cumsum()
dep1990 = dep[dep.year == 1990][["cell_lat", "cell_lon", "totN_kg_ha"]].rename(columns={"totN_kg_ha": "dep1990"})
W = W.merge(dep, on=["cell_lat", "cell_lon", "year"], how="left").merge(dep1990, on=["cell_lat", "cell_lon"], how="left")
W = W[W.year.between(1990, 2025) & W.dep0.notna()]
if SRC == "core":
    for g in ["REF", "SENS", "TOL", "CANTH", "ECM_ALL"]:
        W[g + "_noinat"] = W[g]
W.to_parquet(f"out/{SRC}_nitrogen_panel.parquet")

res = {}
# 1. Deposition history, area-weighted by unit, per country
dc = W.drop_duplicates(["gid", "year"]).groupby(["country", "year"]).apply(
    lambda x: np.average(x.dep0, weights=x.area_km2), include_groups=False).unstack(0)
dc.to_csv(f"out/{SRC}_ndep_country_year.csv")
res["dep_1990"] = dc.loc[1990].round(1).to_dict()
res["dep_2023"] = dc.loc[2023].round(1).to_dict()
res["dep_change_pct"] = ((dc.loc[2023] / dc.loc[1990] - 1) * 100).round(0).to_dict()

# 2. Sensitive share among indicator ECM records, by country and period
def share(df, s="SENS", t="TOL"):
    return df[s].sum() / (df[s].sum() + df[t].sum())

W["period"] = pd.cut(W.year, [1989, 1999, 2009, 2017, 2025], labels=["1990s", "2000s", "2010-17", "2018-25"])
tab = W.groupby(["country", "period"], observed=True).apply(lambda x: share(x), include_groups=False).unstack()
tab2 = W.groupby(["country", "period"], observed=True).apply(lambda x: share(x, "SENS_noinat", "TOL_noinat"), include_groups=False).unstack()
tab.to_csv(f"out/{SRC}_sens_share_country_period.csv")
tab2.to_csv(f"out/{SRC}_sens_share_country_period_noinat.csv")
res["share_by_period"] = tab.round(3).to_dict()
res["share_by_period_noinat"] = tab2.round(3).to_dict()

# 3. Fixed-effects Poisson: sensitive records with exposure = all indicator records
P = W[(W.SENS + W.TOL) >= 10].copy()
P["ind"] = P.SENS + P.TOL
P["log_ind"] = np.log(P.ind)
P["cy"] = P.country + "_" + P.year.astype(str)
models = {}
for v in ["dep0", "dep_ma3", "dep_ma5", "dep_ma10", "dep_lag5", "dep_lag10"]:
    Q = P[P[v].notna()].copy()
    m = pf.fepois(f"SENS ~ {v} | gid + cy", data=Q, offset="log_ind", vcov={"CRV1": "gid"}, fixef_maxiter=200000, fixef_tol=1e-7)
    t = m.tidy()
    models[v] = dict(coef=float(t.loc[v, "Estimate"]), se=float(t.loc[v, "Std. Error"]),
                     p=float(t.loc[v, "Pr(>|t|)"]), n=int(m._N), deviance=float(m.deviance))
res["fe_poisson"] = models

# same model without iNaturalist records
Pn = W[(W.SENS_noinat + W.TOL_noinat) >= 10].copy()
Pn["log_ind"] = np.log(Pn.SENS_noinat + Pn.TOL_noinat)
Pn["cy"] = Pn.country + "_" + Pn.year.astype(str)
m = pf.fepois("SENS_noinat ~ dep_ma5 | gid + cy", data=Pn, offset="log_ind", vcov={"CRV1": "gid"}, fixef_maxiter=200000, fixef_tol=1e-7)
t = m.tidy()
res["fe_poisson_noinat_ma5"] = dict(coef=float(t.loc["dep_ma5", "Estimate"]), se=float(t.loc["dep_ma5", "Std. Error"]),
                                    p=float(t.loc["dep_ma5", "Pr(>|t|)"]), n=int(m._N))

# placebo: saprotroph-dominated reference share should not respond the same way.
# Tolerant genera against the whole Agaricomycetes reference
Pt = W[(W.REF >= 50) & (W.TOL >= 1)].copy()
Pt["log_ref"] = np.log(Pt.REF)
Pt["cy"] = Pt.country + "_" + Pt.year.astype(str)
for resp in ["TOL", "SENS"]:
    m = pf.fepois(f"{resp} ~ dep_ma5 | gid + cy", data=Pt, offset="log_ref", vcov={"CRV1": "gid"}, fixef_maxiter=200000, fixef_tol=1e-7)
    t = m.tidy()
    res[f"fe_poisson_{resp}_vs_REF"] = dict(coef=float(t.loc["dep_ma5", "Estimate"]), se=float(t.loc["dep_ma5", "Std. Error"]),
                                            p=float(t.loc["dep_ma5", "Pr(>|t|)"]), n=int(m._N))

# 4. Long difference: change in logit share 1990-1999 to 2016-2025 against deposition fall and legacy
def logit(p):
    return np.log(p / (1 - p))

a = W[W.year.between(1990, 1999)].groupby("gid")[["SENS", "TOL"]].sum()
b = W[W.year.between(2016, 2025)].groupby("gid")[["SENS", "TOL"]].sum()
LD = a.join(b, lsuffix="_a", rsuffix="_b")
LD = LD[(LD.SENS_a + LD.TOL_a >= 30) & (LD.SENS_b + LD.TOL_b >= 30)]
LD["s_a"] = (LD.SENS_a + 0.5) / (LD.SENS_a + LD.TOL_a + 1)
LD["s_b"] = (LD.SENS_b + 0.5) / (LD.SENS_b + LD.TOL_b + 1)
LD["dlogit"] = logit(LD.s_b) - logit(LD.s_a)
da = W[W.year.between(1990, 1999)].groupby("gid").dep0.mean()
db = W[W.year.between(2014, 2023)].groupby("gid").dep0.mean()
LD = LD.join(da.rename("dep_a")).join(db.rename("dep_b")).join(W.groupby("gid").country.first())
LD["ddep"] = LD.dep_b - LD.dep_a
LD["w"] = 1 / (1 / (LD.SENS_a + 0.5) + 1 / (LD.TOL_a + 0.5) + 1 / (LD.SENS_b + 0.5) + 1 / (LD.TOL_b + 0.5))
LD.to_csv(f"out/{SRC}_nitrogen_long_difference.csv")
import statsmodels.formula.api as smf
m1 = smf.wls("dlogit ~ ddep + C(country)", data=LD, weights=LD.w).fit(cov_type="HC1")
m2 = smf.wls("dlogit ~ dep_a + dep_b + C(country)", data=LD, weights=LD.w).fit(cov_type="HC1")
res["long_diff"] = dict(n=int(len(LD)), ddep_coef=float(m1.params["ddep"]), ddep_se=float(m1.bse["ddep"]), ddep_p=float(m1.pvalues["ddep"]),
                        dep_a_coef=float(m2.params["dep_a"]), dep_a_p=float(m2.pvalues["dep_a"]),
                        dep_b_coef=float(m2.params["dep_b"]), dep_b_p=float(m2.pvalues["dep_b"]),
                        share_units_improved=float((LD.dlogit > 0).mean()),
                        median_dlogit=float(LD.dlogit.median()))

# 5. Threshold: binned relationship of sensitive share with current 5-year deposition
P["bin"] = pd.cut(P.dep_ma5, [0, 3, 5, 7.5, 10, 12.5, 15, 20, 25, 40])
bt = P.groupby("bin", observed=True).apply(lambda x: pd.Series({"share": share(x), "unit_years": len(x), "records": int(x.ind.sum())}), include_groups=False)
bt.to_csv(f"out/{SRC}_sens_share_by_deposition_bin.csv")
res["share_by_dep_bin"] = {str(k): v for k, v in bt.round(3).to_dict("index").items()}

# 6. Chanterelle against all Agaricomycetes, by country and period (with and without iNat)
ct = W.groupby(["country", "period"], observed=True).apply(lambda x: 1000 * x.CANTH.sum() / x.REF.sum(), include_groups=False).unstack()
ctn = W.groupby(["country", "period"], observed=True).apply(lambda x: 1000 * x.CANTH_noinat.sum() / x.REF_noinat.sum(), include_groups=False).unstack()
ct.to_csv(f"out/{SRC}_chanterelle_per_1000.csv"); ctn.to_csv(f"out/{SRC}_chanterelle_per_1000_noinat.csv")
res["chanterelle_per_1000"] = ct.round(2).to_dict()
res["chanterelle_per_1000_noinat"] = ctn.round(2).to_dict()

# 7. Exceedance: share of unit area above 10 and 5.5 kg N/ha/yr
ex = W.drop_duplicates(["gid", "year"])
exc = ex.groupby("year").apply(lambda x: pd.Series({
    "area_gt10": np.average(x.dep0 > 10, weights=x.area_km2),
    "area_gt5_5": np.average(x.dep0 > 5.5, weights=x.area_km2)}), include_groups=False)
exc.to_csv(f"out/{SRC}_exceedance_area.csv")
res["exceedance_1990_2023"] = exc.loc[[1990, 2000, 2010, 2023]].round(3).to_dict("index")
json.dump(res, open(f"out/{SRC}_nitrogen_results.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))

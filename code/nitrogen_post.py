"""Dose curve, robust legacy test, curve-versus-observed recovery and projections.
Appends to out/legacy_results.json (run after nitrogen_legacy.py)."""
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

W = pd.read_parquet("out/core_nitrogen_panel.parquet")
L = json.load(open("out/legacy_results.json"))

# dose curve across all unit-years with >= 10 indicator records
P = W[(W.SENS + W.TOL) >= 10].copy()
P["n"] = P.SENS + P.TOL
P["p"] = P.SENS / P.n
m = smf.glm("p ~ np.log(dep_ma5)", data=P, family=sm.families.Binomial(), var_weights=P.n).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(P.gid)[0]})
a, b = [float(v) for v in m.params]
inv = lambda s: float(np.exp((np.log(s / (1 - s)) - a) / b))
f = lambda d: float(1 / (1 + np.exp(-(a + b * np.log(d)))))
L["dose_curve"] = dict(a=a, b=b, se_b=float(m.bse.iloc[1]), load_for_share={str(s): inv(s) for s in [0.5, 0.75, 0.85, 0.9]},
                       n_unit_years=int(len(P)), n_records=int(P.n.sum()))

# robust legacy test, 2014-2025
late = W[W.year.between(2014, 2025)]
g = late.groupby("gid").agg(SENS=("SENS", "sum"), TOL=("TOL", "sum"), dep_now=("dep_ma5", "mean"),
                            dep1990=("dep1990", "first"), country=("country", "first"))
g = g[(g.SENS + g.TOL) >= 30].copy()
g["n"] = g.SENS + g.TOL
g["p"] = g.SENS / g.n
out = {}
for name, sub in [("all", g), ("stable", g[g.country.isin(["GB", "DK", "SE", "BE", "NO"])])]:
    mm = smf.glm("p ~ bs(dep_now, df=4) + dep1990 + C(country)", data=sub, family=sm.families.Binomial(),
                 var_weights=sub.n).fit(cov_type="HC1")
    out[name] = dict(coef=float(mm.params["dep1990"]), p=float(mm.pvalues["dep1990"]),
                     odds_change_per_kg=float(np.exp(mm.params["dep1990"]) - 1), n_units=int(len(sub)))
L["legacy_glm_robust"] = out

# curve versus observed trajectories
cvo = {}
for k, v in list(L.items()):
    if k.startswith("traj_trend_"):
        c = k[len("traj_trend_"):]
        cvo[c] = dict(dep=[v["dep_first"], v["dep_last"]], obs=[v["share_first5"], v["share_last5"]],
                      pred=[f(v["dep_first"]), f(v["dep_last"])])
L["curve_vs_observed"] = cvo

# projections on fixed-area country deposition
dc = pd.read_csv("out/ndep_country_year_fixed.csv", index_col=0)
t75, t85 = inv(0.75), inv(0.85)
proj = {}
for c in dc.columns:
    s = dc[c].loc[2008:2023]
    sl = np.polyfit(s.index, np.log(s.values), 1)[0]
    now = float(dc[c].loc[2023])
    proj[c] = dict(now=now, pct_per_year=float((np.exp(sl) - 1) * 100), pred_share=f(now),
                   year_75=(int(round(2023 + np.log(t75 / now) / sl)) if now > t75 else None),
                   year_85=(int(round(2023 + np.log(t85 / now) / sl)) if now > t85 else None))
L["projection_fixed"] = proj
json.dump(L, open("out/legacy_results.json", "w"), indent=1, default=str)
print(json.dumps({k: L[k] for k in ["dose_curve", "legacy_glm_robust", "curve_vs_observed", "projection_fixed"]}, indent=1, default=str))

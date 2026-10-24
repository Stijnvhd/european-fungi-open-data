"""Legacy test: at the same current nitrogen load, do places that were once heavily
polluted still hold fewer nitrogen-sensitive ectomycorrhizal fungi than places that
were never heavily polluted? Run on the core recording-scheme panel."""
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

W = pd.read_parquet("out/core_nitrogen_panel.parquet")
STABLE = ["GB", "DK", "SE", "BE", "NO"]  # schemes recording consistently since the 1990s
res = {}

late = W[W.year.between(2014, 2025)].copy()
g = late.groupby("gid").agg(SENS=("SENS", "sum"), TOL=("TOL", "sum"), dep_now=("dep_ma5", "mean"),
                            dep1990=("dep1990", "first"), country=("country", "first"))
g = g[(g.SENS + g.TOL) >= 30].copy()
g["ind"] = g.SENS + g.TOL
g["legacy"] = g.dep1990 - g.dep_now
g["share"] = g.SENS / g.ind
g["now_bin"] = pd.cut(g.dep_now, [0, 5, 7.5, 10, 12.5, 15, 20, 40])
g["was_high"] = np.where(g.dep1990 >= 20, "1990 load 20 kg or more", "1990 load under 20 kg")
tab = g.groupby(["now_bin", "was_high"], observed=True).apply(
    lambda x: pd.Series({"share": x.SENS.sum() / x.ind.sum(), "units": len(x), "records": int(x.ind.sum())}), include_groups=False)
tab.to_csv("out/legacy_table.csv")
print(tab)
res["legacy_table"] = {f"{a}|{b}": v for (a, b), v in tab.round(3).to_dict("index").items()}

# binomial GLM: share ~ current load (spline) + legacy + country
g["log_ind"] = np.log(g.ind)
m = smf.glm("SENS ~ bs(dep_now, df=4) + legacy + C(country)", data=g, family=sm.families.Binomial(),
            var_weights=None, freq_weights=None, exposure=None, offset=None,
            ).fit()
# use counts properly: endog as (SENS, TOL)
y = np.column_stack([g.SENS, g.TOL])
X = smf.glm("SENS ~ bs(dep_now, df=4) + legacy + C(country)", data=g).exog
names = smf.glm("SENS ~ bs(dep_now, df=4) + legacy + C(country)", data=g).exog_names
mb = sm.GLM(y, X, family=sm.families.Binomial()).fit(scale="X2")
i = names.index("legacy")
res["legacy_glm"] = dict(coef_per_kg=float(mb.params[i]), se=float(mb.bse[i]), p=float(mb.pvalues[i]),
                         odds_ratio_per_10kg=float(np.exp(10 * mb.params[i])), n_units=int(len(g)), dispersion=float(mb.scale))
print(res["legacy_glm"])

# yearly trajectories by 1990 deposition class, stable schemes only
S = W[W.country.isin(STABLE)].copy()
S["cls"] = pd.cut(S.dep1990, [0, 10, 20, 60], labels=["under 10 kg in 1990", "10 to 20 kg", "over 20 kg"])
traj = S.groupby(["cls", "year"], observed=True).apply(
    lambda x: pd.Series({"share": x.SENS.sum() / (x.SENS.sum() + x.TOL.sum()), "ind": int(x.SENS.sum() + x.TOL.sum()),
                         "dep": float(np.average(x.dep0, weights=x.area_km2))}), include_groups=False).reset_index()
traj.to_csv("out/legacy_trajectories.csv", index=False)
for c in traj.cls.unique():
    t = traj[(traj.cls == c) & (traj.ind >= 100)]
    mm = smf.wls("np.log(share/(1-share)) ~ I(year-2000)", data=t, weights=np.sqrt(t.ind)).fit()
    res[f"traj_trend_{c}"] = dict(logit_per_decade=float(mm.params.iloc[1] * 10), p=float(mm.pvalues.iloc[1]),
                                  share_first5=float(t.share.iloc[:5].mean()), share_last5=float(t.share.iloc[-5:].mean()),
                                  dep_first=float(t.dep.iloc[:5].mean()), dep_last=float(t.dep.iloc[-5:].mean()),
                                  years=f"{t.year.min()}-{t.year.max()}")
print({k: v for k, v in res.items() if k.startswith("traj")})

# stable-scheme long difference by 1990 class
a = S[S.year.between(1990, 1999)].groupby("gid")[["SENS", "TOL"]].sum()
b = S[S.year.between(2016, 2025)].groupby("gid")[["SENS", "TOL"]].sum()
LD = a.join(b, lsuffix="_a", rsuffix="_b").join(S.groupby("gid").cls.first())
LD = LD[(LD.SENS_a + LD.TOL_a >= 30) & (LD.SENS_b + LD.TOL_b >= 30)]
LD["s_a"] = LD.SENS_a / (LD.SENS_a + LD.TOL_a)
LD["s_b"] = LD.SENS_b / (LD.SENS_b + LD.TOL_b)
res["stable_ld_by_class"] = LD.groupby("cls", observed=True).apply(
    lambda x: {"units": int(len(x)), "share_1990s": float(x.SENS_a.sum() / (x.SENS_a.sum() + x.TOL_a.sum())),
               "share_2016_25": float(x.SENS_b.sum() / (x.SENS_b.sum() + x.TOL_b.sum())),
               "units_improved": float((x.s_b > x.s_a).mean())}, include_groups=False).to_dict()
print(res["stable_ld_by_class"])
json.dump(res, open("out/legacy_results.json", "w"), indent=1, default=str)

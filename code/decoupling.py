"""Are the fungi that live with trees keeping time with their trees?
Fungal season dates (2001-2023) against MODIS land surface phenology at the same sites."""
import json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

MO = pd.read_parquet("data/site_modis.parquet")
MO = MO[MO.npix >= 10]
L = MO.groupby(["country", "year", "band"]).doy_median.median().unstack("band").reset_index()
L = L.rename(columns={"MidGreenup": "leaf_up", "MidGreendown": "leaf_down", "Dormancy": "leaf_dorm"})
L.to_csv("out/country_leaf_phenology.csv", index=False)
K = pd.read_csv("out/country_climate.csv")
M = pd.read_parquet("out/phenology_metrics.parquet")
res = {}

# 1. trends in leaf dates 2001-2023
for v in ["leaf_up", "leaf_down", "leaf_dorm"]:
    d = L.dropna(subset=[v])
    m = smf.ols(f"{v} ~ I((year-2010)/10) + C(country)", data=d).fit(cov_type="HC1")
    res[f"trend_{v}"] = dict(per_decade=float(m.params.iloc[-1]), se=float(m.bse.iloc[-1]), p=float(m.pvalues.iloc[-1]),
                             mean=float(d[v].mean()))
    pc = {}
    for c, g in d.groupby("country"):
        mm = smf.ols(f"{v} ~ I((year-2010)/10)", data=g).fit()
        pc[c] = dict(per_decade=float(mm.params.iloc[1]), p=float(mm.pvalues.iloc[1]), mean=float(g[v].mean()))
    res[f"trend_{v}_by_country"] = pc

# 2. fungal trends over the same window, and the gap
for ver in ["effort", "core"]:
    d = M[(M.version == ver) & M.year.between(2001, 2023)].merge(L, on=["country", "year"]).merge(K, on=["country", "year"])
    d["sc"] = d.taxon + "_" + d.country
    d["w"] = np.sqrt(d.n.clip(upper=2000))
    d["yr"] = (d.year - 2010) / 10
    d["gap_end"] = d["end"] - d.leaf_down
    d["gap_mid"] = d["middle"] - d.leaf_down
    for gu in ["ECM", "SAP"]:
        dd = d[d.guild == gu]
        out = {}
        for var in ["end", "middle", "gap_end", "gap_mid"]:
            m = smf.wls(f"{var} ~ yr + C(sc)", data=dd, weights=dd.w).fit(cov_type="cluster", cov_kwds={"groups": dd.sc})
            out[f"trend_{var}"] = dict(per_decade=float(m.params["yr"]), se=float(m.bse["yr"]), p=float(m.pvalues["yr"]))
        out["mean_gap_end"] = float(np.average(dd.gap_end, weights=dd.w))
        out["mean_gap_mid"] = float(np.average(dd.gap_mid, weights=dd.w))
        # 3. interannual coupling: fungal date on leaf-down date, with species-country effects and a common year trend,
        #    with and without autumn temperature
        for var in ["end", "middle"]:
            m1 = smf.wls(f"{var} ~ leaf_down + yr + C(sc)", data=dd, weights=dd.w).fit(cov_type="cluster", cov_kwds={"groups": dd.sc})
            m2 = smf.wls(f"{var} ~ leaf_down + t_so + t_on + t_jja + p_as + yr + C(sc)", data=dd, weights=dd.w).fit(
                cov_type="cluster", cov_kwds={"groups": dd.sc})
            out[f"coupling_{var}"] = dict(b=float(m1.params["leaf_down"]), se=float(m1.bse["leaf_down"]), p=float(m1.pvalues["leaf_down"]),
                                          b_ctrl=float(m2.params["leaf_down"]), se_ctrl=float(m2.bse["leaf_down"]),
                                          p_ctrl=float(m2.pvalues["leaf_down"]), n=int(len(dd)))
        res[f"{ver}|{gu}"] = out
    # guild difference in coupling, interaction model
    d["ecm"] = (d.guild == "ECM").astype(int)
    m = smf.wls("end ~ leaf_down*ecm + t_so + t_on + t_jja + p_as + yr + C(sc)", data=d, weights=d.w).fit(
        cov_type="cluster", cov_kwds={"groups": d.sc})
    res[f"{ver}|coupling_end_interaction"] = dict(b=float(m.params["leaf_down:ecm"]), se=float(m.bse["leaf_down:ecm"]), p=float(m.pvalues["leaf_down:ecm"]))
    # guild difference in trend of gap
    m = smf.wls("gap_end ~ yr*ecm + C(sc)", data=d, weights=d.w).fit(cov_type="cluster", cov_kwds={"groups": d.sc})
    res[f"{ver}|gap_trend_interaction"] = dict(b=float(m.params["yr:ecm"]), se=float(m.bse["yr:ecm"]), p=float(m.pvalues["yr:ecm"]))
json.dump(res, open("out/decoupling_results.json", "w"), indent=1)
print(json.dumps(res, indent=1))

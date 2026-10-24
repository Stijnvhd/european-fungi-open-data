"""The fungal clock: autumn fruiting season of 50 common species in nine European
countries, 1980-2025, from dekad-resolution GBIF counts.

Metrics per species x country x year (records between 1 July and 31 December):
  start = 10th percentile date, middle = median date, end = 90th percentile date,
  length = end - start. Percentiles come from the cumulative record curve, with
  records spread evenly within each dekad.

Three versions of every metric:
  raw        all records
  effort     each dekad's records divided by that dekad's share of all
             vascular plant records in the same country and year (observer
             activity index), so a shift in when people go out does not
             masquerade as a shift in when fungi appear
  core       records of each country's main long-running recording scheme only
"""
import json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

A = pd.read_parquet("data/dekads_fb.parquet")
CORE = pd.read_parquet("data/dekads_core.parquet") if __import__("os").path.exists("data/dekads_core.parquet") else None

def dekad_edges(m, i):
    start = pd.Timestamp(2001, m, 1 + 10 * i).dayofyear
    if i < 2:
        end = start + 10
    else:
        end = (pd.Timestamp(2001, m, 1) + pd.offsets.MonthEnd(0)).dayofyear + 1
    return start, end

ED = {(m, i): dekad_edges(m, i) for m in range(1, 13) for i in range(3)}
for df in [A] + ([CORE] if CORE is not None else []):
    df["d0"] = [ED[(m, i)][0] for m, i in zip(df.month, df.dekad)]
    df["d1"] = [ED[(m, i)][1] for m, i in zip(df.month, df.dekad)]

# observer activity index from vascular plant records
plants = A[A.taxon == "Tracheophyta"][["country", "year", "month", "dekad", "n"]].rename(columns={"n": "eff"})
tot = plants.groupby(["country", "year"]).eff.transform("sum")
plants["effshare"] = plants.eff / tot * 36.0  # 1 = average dekad
A = A.merge(plants[["country", "year", "month", "dekad", "effshare"]], on=["country", "year", "month", "dekad"], how="left")

def quantiles(g, w):
    g = g.sort_values("d0")
    c = np.cumsum(w) / w.sum()
    out = {}
    for q in [0.1, 0.5, 0.9]:
        k = int(np.searchsorted(c.values, q))
        k = min(k, len(g) - 1)
        prev = c.values[k - 1] if k > 0 else 0.0
        frac = (q - prev) / max(c.values[k] - prev, 1e-12)
        out[q] = g.d0.values[k] + frac * (g.d1.values[k] - g.d0.values[k])
    mean = np.average((g.d0 + g.d1) / 2, weights=w)
    return out[0.1], out[0.5], out[0.9], mean

def metrics(df, version):
    rows = []
    sub = df[(df.guild.isin(["ECM", "SAP"])) & (df.month >= 7) & (df.year <= 2025)]
    for (t, gu, c, y), g in sub.groupby(["taxon", "guild", "country", "year"]):
        if version == "effort":
            es = g.effshare.where(g.effshare > 0.05)
            if es.isna().any():
                continue
            w = g.n / es
        else:
            w = g.n.astype(float)
        if g.n.sum() < 30:
            continue
        q10, q50, q90, mean = quantiles(g, w)
        rows.append(dict(taxon=t, guild=gu, country=c, year=y, n=int(g.n.sum()), start=q10, middle=q50, end=q90, mean=mean, length=q90 - q10))
    out = pd.DataFrame(rows)
    out["version"] = version
    return out

M = [metrics(A, "raw"), metrics(A, "effort")]
if CORE is not None:
    M.append(metrics(CORE, "core"))
M = pd.concat(M)
M.to_parquet("out/phenology_metrics.parquet")

res = {"n_species_country_years": M.groupby("version").size().to_dict()}

def trend(df, var, by_guild=True, y0=1981):
    d = df[df.year >= y0].copy()
    d["sc"] = d.taxon + "_" + d.country
    d["w"] = np.sqrt(d.n.clip(upper=2000))
    d["yr"] = (d.year - 2000) / 10.0  # per decade
    out = {}
    for gu in (["ECM", "SAP"] if by_guild else ["ALL"]):
        dd = d if gu == "ALL" else d[d.guild == gu]
        if len(dd) < 50:
            continue
        m = smf.wls(f"{var} ~ yr + C(sc)", data=dd, weights=dd.w).fit(cov_type="cluster", cov_kwds={"groups": dd.sc})
        out[gu] = dict(per_decade=float(m.params["yr"]), se=float(m.bse["yr"]), p=float(m.pvalues["yr"]), n=int(len(dd)))
    return out

for v in ["raw", "effort", "core"]:
    if v not in M.version.unique():
        continue
    mv = M[M.version == v]
    res[f"trend_{v}"] = {var: trend(mv, var) for var in ["start", "middle", "end", "length", "mean"]}
    res[f"trend_{v}_since1995"] = {var: trend(mv, var, y0=1995) for var in ["start", "end", "length"]}

# per-country trends (effort version)
me = M[M.version == "effort"]
pc = {}
for c in sorted(me.country.unique()):
    pc[c] = {var: trend(me[me.country == c], var) for var in ["start", "end", "length"]}
res["trend_effort_by_country"] = pc

# per-species slopes (effort, pooled across countries)
sp = []
for (t, gu), g in me[me.year >= 1981].groupby(["taxon", "guild"]):
    if g.year.nunique() < 15:
        continue
    g = g.copy(); g["yr"] = (g.year - 2000) / 10; g["w"] = np.sqrt(g.n.clip(upper=2000))
    for var in ["start", "end", "length"]:
        m = smf.wls(f"{var} ~ yr + C(country)", data=g, weights=g.w).fit()
        sp.append(dict(taxon=t, guild=gu, var=var, slope=m.params["yr"], p=m.pvalues["yr"], n=len(g)))
SP = pd.DataFrame(sp)
SP.to_csv("out/phenology_species_slopes.csv", index=False)
res["species_slopes_summary"] = SP.groupby(["var", "guild"]).apply(
    lambda x: {"n_species": int(len(x)), "share_positive": float((x.slope > 0).mean()),
               "share_sig_pos": float(((x.slope > 0) & (x.p < 0.05)).mean()),
               "share_sig_neg": float(((x.slope < 0) & (x.p < 0.05)).mean()),
               "median_slope": float(x.slope.median())}, include_groups=False)
res["species_slopes_summary"] = {f"{a}|{b}": v for (a, b), v in res["species_slopes_summary"].items()}

# observer seasonality itself: the plant-record centre of gravity Jul-Dec over time
pl = A[(A.taxon == "Tracheophyta") & (A.month >= 7) & (A.year <= 2025)]
pcg = pl.groupby(["country", "year"]).apply(lambda g: np.average((g.d0 + g.d1) / 2, weights=g.n) if g.n.sum() > 0 else np.nan, include_groups=False).rename("plant_centre").reset_index()
pcg.to_csv("out/plant_record_centre.csv", index=False)
m = smf.ols("plant_centre ~ I((year-2000)/10) + C(country)", data=pcg[pcg.year >= 1981]).fit(cov_type="HC1")
res["observer_autumn_centre_shift_per_decade"] = dict(coef=float(m.params.iloc[-1]), p=float(m.pvalues.iloc[-1]))

# spring fruiting share (Mar-Jun) among records, by guild and decade
sub = A[A.guild.isin(["ECM", "SAP"]) & (A.year <= 2025)].copy()
sub["spring"] = sub.month.between(3, 6)
sub["decade"] = (sub.year // 10) * 10
sp_share = sub.groupby(["guild", "decade"]).apply(lambda g: g.n[g.spring].sum() / g.n.sum(), include_groups=False).unstack(0)
sp_share.to_csv("out/spring_share.csv")
res["spring_share_by_decade"] = {g: {str(k): v for k, v in d.items()} for g, d in sp_share.round(4).to_dict().items()}

json.dump(res, open("out/phenology_results.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))

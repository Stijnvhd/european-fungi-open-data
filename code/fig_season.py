import sys
sys.path.insert(0, "scripts")
from figstyle import *
import json
import numpy as np
import pandas as pd

M = pd.read_parquet("out/phenology_metrics.parquet")
R = json.load(open("out/phenology_results.json"))

def anomalies(ver):
    d = M[(M.version == ver) & (M.year.between(1985, 2025))].copy()
    d["sc"] = d.taxon + "_" + d.country
    d["w"] = np.sqrt(d.n.clip(upper=2000))
    for v in ["start", "middle", "end"]:
        d[v + "_a"] = d[v] - d.groupby("sc")[v].transform("mean")
    g = d.groupby(["guild", "year"]).apply(
        lambda x: pd.Series({v: np.average(x[v + "_a"], weights=x.w) for v in ["start", "middle", "end"]}), include_groups=False)
    return g.reset_index()

A = anomalies("effort")
A.to_csv("out/season_anomalies_effort.csv", index=False)
fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.5), sharey=True)
for ax, gu, title in [(axs[0], "ECM", "Fungi that live with trees"), (axs[1], "SAP", "Fungi that rot wood and leaves")]:
    a = A[A.guild == gu].set_index("year")
    for v, col, lab in [("end", ORANGE, "end (90% of records)"), ("middle", BLUE, "middle (median)"), ("start", AQUA, "start (10% of records)")]:
        s = a[v].rolling(5, center=True, min_periods=3).mean()
        ax.plot(a.index, a[v], color=col, lw=0.8, alpha=0.35)
        ax.plot(s.index, s.values, color=col, lw=2.2, label=lab)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_title(title); ax.set_xlim(1985, 2025); ax.set_xticks([1990, 2000, 2010, 2020])
axs[0].set_ylabel("days later (+) or earlier (-) than average")
axs[1].legend(loc="upper left", fontsize=7.2)
fig.subplots_adjust(wspace=0.08, bottom=0.17)
fig.text(0.07, 0.02, "Autumn season dates of 50 species in nine countries, relative to each species and country's average; effort-corrected. Thin: single years; thick: five-year mean.", fontsize=6.6, color=MUTED)
save(fig, "season_anomalies")

# per-country trends: start and end, both guilds pooled by averaging the two guild estimates
names = {"NL": "Netherlands", "BE": "Belgium", "DE": "Germany", "DK": "Denmark", "FR": "France", "GB": "Britain", "SE": "Sweden", "FI": "Finland", "NO": "Norway"}
rows = []
for c, d in R["trend_effort_by_country"].items():
    for var in ["start", "end"]:
        for gu, x in d[var].items():
            rows.append(dict(country=c, var=var, guild=gu, b=x["per_decade"], se=x["se"]))
T = pd.DataFrame(rows)
T.to_csv("out/season_trend_by_country.csv", index=False)
order = T[(T["var"] == "end")].groupby("country").b.mean().sort_values().index
fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.6), sharey=True)
for ax, var, title in [(axs[0], "start", "Start of season"), (axs[1], "end", "End of season")]:
    for k, (gu, col, off) in enumerate([("ECM", BLUE, -0.15), ("SAP", ORANGE, 0.15)]):
        t = T[(T["var"] == var) & (T.guild == gu)].set_index("country").reindex(order)
        y = np.arange(len(order)) + off
        ax.errorbar(t.b, y, xerr=1.96 * t.se, fmt="o", color=col, ms=5, lw=1.2, capsize=0,
                    label="live with trees" if gu == "ECM" else "rot wood and leaves")
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(np.arange(len(order))); ax.set_yticklabels([names[c] for c in order]); ax.grid(axis="y", visible=False)
    ax.set_title(title); ax.set_xlabel("days per decade")
axs[1].legend(loc="lower right", fontsize=7.2)
fig.subplots_adjust(wspace=0.08, bottom=0.2)
fig.text(0.18, 0.02, "Trend 1981-2025 with 95% interval, effort-corrected. Finland: too few rot-fungus records.", fontsize=7, color=MUTED)
save(fig, "season_by_country")
print(T.pivot_table(index="country", columns=["var", "guild"], values="b").round(1))

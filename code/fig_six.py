import sys
sys.path.insert(0, "scripts")
from figstyle import *
import json
import numpy as np
import pandas as pd

C = pd.read_parquet("data/dekads_core.parquet")
C = C[C.year <= 2025]
SIX = [("Cantharellus cibarius", "Chanterelle"), ("Amanita muscaria", "Fly agaric"), ("Boletus edulis", "Cep"),
       ("Hypholoma fasciculare", "Sulphur tuft"), ("Armillaria mellea", "Honey fungus"), ("Flammulina velutipes", "Velvet shank")]
P1, P2 = (1986, 2000), (2011, 2025)
C["k"] = (C.month - 1) * 3 + C.dekad  # 0..35
out = {}
fig, axs = plt.subplots(3, 2, figsize=(6.8, 6.6), sharex=True)
for ax, (sp, lab) in zip(axs.ravel(), SIX):
    d = C[C.taxon == sp]
    curves = {}
    for (a, b), col, name in [(P1, "#9fb9d9", f"{P1[0]}-{P1[1]}"), (P2, ORANGE, f"{P2[0]}-{P2[1]}")]:
        x = d[d.year.between(a, b)]
        # weight countries equally: normalise each country's distribution, then average
        per = x.groupby(["country", "k"]).n.sum().unstack(fill_value=0).reindex(columns=range(36), fill_value=0)
        per = per[per.sum(axis=1) >= 50]
        dist = per.div(per.sum(axis=1), axis=0).mean()
        curves[name] = dist
        ax.fill_between(range(36), dist.values * 100, color=col, alpha=0.35 if col != ORANGE else 0.0, lw=0)
        ax.plot(range(36), dist.values * 100, color=col if col == ORANGE else BLUE, lw=1.8, label=name)
    # summary: median and 90th percentile day, autumn window
    summ = {}
    for name, dist in curves.items():
        aut = dist.iloc[18:]
        cdf = (aut / aut.sum()).cumsum().values
        mids = np.array([(pd.Timestamp(2001, m, 1) + pd.Timedelta(days=4 + 10 * i)).dayofyear for m in range(7, 13) for i in range(3)])
        summ[name] = dict(median=float(np.interp(0.5, cdf, mids)), p90=float(np.interp(0.9, cdf, mids)),
                          p10=float(np.interp(0.1, cdf, mids)), countries=int(len(d.country.unique())))
    out[sp] = summ
    ax.set_title(lab, fontsize=9.5)
    ax.set_xticks([0, 6, 12, 18, 24, 30]); ax.set_xticklabels(["Jan", "Mar", "May", "Jul", "Sep", "Nov"])
    ax.set_yticks([])
axs[0, 0].legend(loc="upper left", fontsize=7.2)
fig.subplots_adjust(hspace=0.45, wspace=0.12, bottom=0.08)
fig.text(0.07, 0.015, "Share of each species' records in each ten-day period of the year, long-running national schemes, each country weighted equally.", fontsize=7, color=MUTED)
save(fig, "six_seasons")
json.dump(out, open("out/six_seasons.json", "w"), indent=1)
for k, v in out.items():
    print(k, {n: {a: round(b) for a, b in s.items()} for n, s in v.items()})

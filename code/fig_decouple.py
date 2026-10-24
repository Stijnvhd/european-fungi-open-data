import sys
sys.path.insert(0, "scripts")
from figstyle import *
import numpy as np
import pandas as pd

L = pd.read_csv("out/country_leaf_phenology.csv")
A = pd.read_csv("out/season_anomalies_effort.csv")
la = L.copy()
for v in ["leaf_up", "leaf_down"]:
    la[v] = la[v] - la.groupby("country")[v].transform("mean")
leaf = la.groupby("year")[["leaf_up", "leaf_down"]].mean()
fung = A[A.year.between(2001, 2023)].pivot(index="year", columns="guild", values="end")
fung = fung - fung.mean()
fig, ax = plt.subplots(figsize=(6.6, 4.0))
for s, col, lab, lw in [(leaf.leaf_up, AQUA, "leaves: spring green-up", 1.6), (leaf.leaf_down, "#8a8985", "leaves: autumn brown-down", 1.6),
                        (fung["ECM"], BLUE, "fungi with trees: end of season", 2.2), (fung["SAP"], ORANGE, "fungi of wood and leaves: end of season", 2.2)]:
    sm = s.rolling(5, center=True, min_periods=3).mean()
    ax.plot(s.index, s.values, color=col, lw=0.7, alpha=0.35)
    ax.plot(sm.index, sm.values, color=col, lw=lw, label=lab)
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlim(2001, 2023); ax.set_ylabel("days later (+) or earlier (-) than 2001-2023 average")
ax.legend(loc="lower left", fontsize=7.2, ncol=2)
ax.set_xticks([2005, 2010, 2015, 2020]); ax.set_ylim(-36, 13)
ax.set_title("The fungi move later; the leaves do not")
source(ax, "Leaf dates: MODIS MCD12Q2 at 108 sites in nine countries. Fungal dates: effort-corrected season end, 50 species. Thin: single years; thick: five-year mean.")
save(fig, "fungi_vs_leaves")

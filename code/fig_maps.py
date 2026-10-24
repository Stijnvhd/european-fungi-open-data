import sys
sys.path.insert(0, "scripts")
from figstyle import *
import glob
import json
import numpy as np
import pandas as pd
from shapely.geometry import shape
from shapely.ops import unary_union
from matplotlib.colors import BoundaryNorm, ListedColormap

# country outlines from GADM level 2, simplified
outl = {}
for f in glob.glob("data/gadm/gadm41_*_2.json"):
    j = json.load(open(f))
    g = unary_union([shape(ft["geometry"]).buffer(0) for ft in j["features"]]).simplify(0.03)
    outl[f.split("_")[1]] = g

def draw_outlines(ax):
    for g in outl.values():
        geoms = getattr(g, "geoms", [g])
        for p in geoms:
            if p.area < 0.02:
                continue
            x, y = p.exterior.xy
            ax.plot(x, y, color="#55534e", lw=0.45)

D = pd.read_parquet("data/emep_ndep_halfdeg.parquet")
lat = np.arange(44, 72, 0.5); lon = np.arange(-11, 32, 0.5)
bounds = [0, 2.5, 5, 7.5, 10, 15, 20, 30, 60]
cols = ["#f3f6fb", "#d6e6f8", "#a9cbf0", "#6fa7e4", "#2a78d6", "#1a56a1", "#0f3768", "#08203f"]
cmap = ListedColormap(cols); norm = BoundaryNorm(bounds, cmap.N)
fig, axs = plt.subplots(1, 2, figsize=(6.8, 4.4))
for ax, y in zip(axs, [1990, 2023]):
    d = D[D.year == y].pivot(index="cell_lat", columns="cell_lon", values="totN_kg_ha").reindex(index=lat, columns=lon)
    m = ax.pcolormesh(lon, lat, d.values, cmap=cmap, norm=norm, shading="auto")
    draw_outlines(ax)
    ax.set_xlim(-8, 31); ax.set_ylim(46, 71.5); ax.set_aspect(1.7)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(str(y))
cb = fig.colorbar(m, ax=axs, orientation="horizontal", fraction=0.05, pad=0.04, ticks=bounds)
cb.set_label("total nitrogen deposition, kg per hectare per year", color=INK2, fontsize=8)
cb.outline.set_visible(False)
fig.text(0.07, 0.01, "EMEP MSC-W rv5.6 trend runs, grid average, aggregated to half-degree cells.", fontsize=7, color=MUTED)
save(fig, "ndep_maps")

# district change in sensitive share, long-running schemes
W = pd.read_parquet("out/core_nitrogen_panel.parquet")
S = W[W.country.isin(["GB", "DK", "SE", "BE", "NO"])]
a = S[S.year.between(1990, 1999)].groupby("gid")[["SENS", "TOL"]].sum()
b = S[S.year.between(2016, 2025)].groupby("gid")[["SENS", "TOL"]].sum()
L = a.join(b, lsuffix="_a", rsuffix="_b")
L = L[(L.SENS_a + L.TOL_a >= 30) & (L.SENS_b + L.TOL_b >= 30)]
L["d"] = (L.SENS_b / (L.SENS_b + L.TOL_b) - L.SENS_a / (L.SENS_a + L.TOL_a)) * 100
cen = pd.read_parquet("data/gadm2_centroids.parquet").set_index("gid")
L = L.join(cen[["lat", "lon"]])
fig, ax = plt.subplots(figsize=(6.6, 5.0))
draw_outlines(ax)
pos = L[L.d >= 0]; neg = L[L.d < 0]
ax.scatter(neg.lon, neg.lat, s=8 + neg.d.abs() * 1.2, color=ORANGE, alpha=0.8, edgecolor="white", linewidth=0.5, label="sensitive share fell")
ax.scatter(pos.lon, pos.lat, s=8 + pos.d.abs() * 1.2, color=BLUE, alpha=0.8, edgecolor="white", linewidth=0.5, label="sensitive share rose")
ax.set_xlim(-8.5, 25); ax.set_ylim(49.3, 69.5); ax.set_aspect(1.7)
ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(loc="lower right", fontsize=7.5)
ax.set_title("Change in the nitrogen-sensitive share, 1990s to 2016-2025")
fig.text(0.07, 0.03, f"{len(L)} districts with at least 30 indicator records in both periods; dot size grows with the size of the change. Long-running schemes of Britain, Denmark, Sweden, Belgium and Norway.", fontsize=6.6, color=MUTED)
save(fig, "district_change_map")
print(len(L), (L.d > 0).mean(), L.d.median())

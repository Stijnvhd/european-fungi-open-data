import sys; sys.path.insert(0,"scripts")
from figstyle import *
import numpy as np, json
S=json.load(open("out/naming_scenarios.json")); r=S["rate_2013_2022"]; n=S["named"]
yrs=np.arange(2023,3100)
fig,ax=plt.subplots(figsize=(6.6,3.5))
for g,col,lab in [(0,MUTED,"rate frozen at 2013-2022 level"),(0.01,AQUA,"rate grows 1% a year"),(0.02,ORANGE,"rate grows 2% a year"),(0.03,BLUE,"rate grows 3% a year")]:
    k=yrs-2022
    cum = n + (r*k if g==0 else r*((1+g)**k-1)/g*(1+g))
    cum=np.where(cum>6.9e6,np.nan,cum); ax.plot(yrs,cum/1e6,color=col,lw=2,label=lab)
for est,lab in [(2.5e6,"best estimate, 2.5 million"),(6.28e6,"molecular estimate, 6.3 million")]:
    ax.axhline(est/1e6,color=INK2,lw=0.8,ls=(0,(3,2))); ax.text(3095,est/1e6+0.1,lab,ha="right",fontsize=7.5,color=INK2)
ax.set_xlim(2023,3100); ax.set_ylim(0,7); ax.set_ylabel("named species (millions)")
ax.legend(loc="center right",bbox_to_anchor=(1,0.55))
ax.set_title("How long to name the fungi")
source(ax,f"Projection from {n:,} accepted species and a 2013-2022 mean of {r:,.0f} new species a year (Species Fungorum Plus, Apr 2024).")
save(fig,"naming_clock")

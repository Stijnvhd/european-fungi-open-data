import sys; sys.path.insert(0,"scripts")
from figstyle import *
import pandas as pd, json
d=pd.read_csv("out/agarico_by_dataset_year.csv")
fig,axs=plt.subplots(1,2,figsize=(6.8,3.3),sharey=False)
spec={"NL":[("8a863029-f435-446a-821e-275f4f641165","Observation.org",ORANGE),("50c9509d-22c7-4a22-a47d-8c48425ef4a7","iNaturalist",YELLOW)],
      "DK":[("84d26682-f762-11e1-a439-00145eb45e9a","Danish Mycological Society",BLUE),("50c9509d-22c7-4a22-a47d-8c48425ef4a7","iNaturalist",YELLOW)]}
for ax,c,title in [(axs[0],"NL","Netherlands"),(axs[1],"DK","Denmark")]:
    x=d[(d.country==c)&(d.year<=2024)]; tot=x[x.dataset=="TOTAL"].set_index("year").n
    parts=[];labs=[];cols=[]
    for k,l,col in spec[c]:
        s=x[x.dataset==k].set_index("year").n.reindex(tot.index,fill_value=0); parts.append(s/1000); labs.append(l); cols.append(col)
    other=(tot-sum(p*1000 for p in parts)).clip(lower=0)/1000
    ax.stackplot(tot.index,*parts,other,colors=cols+["#d9d8d3"],labels=labs+["all other sources"],edgecolor="white",linewidth=0.6)
    ax.set_title(title); ax.set_xlim(1980,2024); ax.set_xticks([1980,2000,2020])
    ax.legend(loc="upper left",fontsize=7.2)
axs[0].set_ylabel("Agaricomycetes records per year (thousands)")
fig.subplots_adjust(wspace=0.28,bottom=0.17)
fig.text(0.07,0.02,"Fruit-body records (observations and specimens) of mushrooms and their relatives on GBIF, by publishing dataset, 1980 to 2024.",fontsize=7,color=MUTED)
save(fig,"record_sources")

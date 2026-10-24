import sys; sys.path.insert(0,"scripts")
from figstyle import *
import pandas as pd, numpy as np, json
a=pd.read_csv("out/descriptions_per_year.csv",index_col=0).iloc[:,0]
a=a.loc[1753:2022]
r=a.rolling(5,center=True,min_periods=3).mean()
fig,ax=plt.subplots(figsize=(6.6,3.4))
ax.bar(a.index,a.values,width=1.0,color="#c9dcf3",linewidth=0)
ax.plot(r.index,r.values,color=BLUE)
for x,lab,y in [(1821,"Fries, Systema\nMycologicum",900),(1882,"Saccardo, Sylloge\nFungorum begins",1500),(1945,"Second World War",300),(2013,"Name registration\nbecomes mandatory",2950)]:
    ax.axvline(x,color=MUTED,lw=0.7,ls=(0,(2,2)))
    ax.text(x+2,y,lab,fontsize=7.2,color=INK2,va="bottom")
ax.set_xlim(1753,2023); ax.set_ylabel("species first named per year")
ax.set_title("Two and a half centuries of naming fungi")
source(ax,"Earliest name among each accepted species and its synonyms. Bars: single years; line: five-year mean. Data: Species Fungorum Plus (Apr 2024), 156,028 species.")
save(fig,"naming_per_year")

T=pd.read_parquet("data/sf_species_years.parquet")
c=T.groupby(T.first_year.astype(int)).size().reindex(range(1753,2024),fill_value=0).cumsum()
fig,ax=plt.subplots(figsize=(6.6,3.4))
ax.fill_between(c.index,0,c.values/1e6,color=BLUE,alpha=0.9,lw=0)
for est,lab in [(1.5e6,"1.5 million (Hawksworth 1991)"),(2.5e6,"2.5 million (Niskanen et al. 2023)"),(6.28e6,"6.3 million molecular taxa (Baldrian et al. 2022)")]:
    ax.axhline(est/1e6,color=ORANGE if est==2.5e6 else MUTED,lw=1.2 if est==2.5e6 else 0.8,ls="-" if est==2.5e6 else (0,(3,2)))
    ax.text(1758,est/1e6+0.08,lab,fontsize=7.5,color=INK2)
ax.text(2022,c.iloc[-1]/1e6+0.12,f"named so far:\n{c.iloc[-1]:,}",ha="right",fontsize=7.5,color=INK)
ax.set_ylim(0,6.8); ax.set_xlim(1753,2024); ax.set_ylabel("species (millions)")
ax.set_title("The named fungi against three estimates of the total")
source(ax,"Cumulative accepted species by year first named. Data: Species Fungorum Plus (Apr 2024); estimates from the cited studies.")
save(fig,"named_vs_estimates")

U=json.load(open("out/unite_summary.json"))
bc=pd.DataFrame(U["by_class"]); bc=bc[bc["class"]!="(unassigned)"].head(12).sort_values("n")
fig,ax=plt.subplots(figsize=(6.6,3.8))
y=np.arange(len(bc))
ax.barh(y,bc.n,color="#d9d8d3",height=0.7,label="species hypotheses without a species name")
ax.barh(y,bc.n*bc.named,color=BLUE,height=0.7,label="with a species name")
ax.set_yticks(y); ax.set_yticklabels(bc["class"]); ax.grid(axis="y",visible=False)
for yi,(n,s) in enumerate(zip(bc.n,bc.named)): ax.text(n+300,yi,f"{s:.0%} named",va="center",fontsize=7.5,color=INK2)
ax.set_xlabel("DNA species hypotheses (UNITE 10, dynamic threshold)")
ax.legend(loc="lower right")
ax.set_title("Most species seen in DNA have no name")
source(ax,f"{U['n_SH']:,} fungal species hypotheses; {U['share_with_species']:.0%} carry a Linnaean species name. Data: UNITE general release 04.04.2024 via SBDI.")
save(fig,"unite_named_by_class")

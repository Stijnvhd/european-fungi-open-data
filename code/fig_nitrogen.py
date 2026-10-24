import sys; sys.path.insert(0,"scripts")
from figstyle import *
import pandas as pd, numpy as np, json
dc=pd.read_csv("out/ndep_country_year_fixed.csv",index_col=0).loc[1990:2023]
names={"NL":"Netherlands","BE":"Belgium","DE":"Germany","DK":"Denmark","FR":"France","GB":"Britain","SE":"Sweden","FI":"Finland","NO":"Norway"}
fig,ax=plt.subplots(figsize=(6.6,3.6))
hi={"NL":ORANGE,"BE":BLUE}
for c in dc.columns:
    col=hi.get(c,"#b9b8b3"); ax.plot(dc.index,dc[c],color=col,lw=2.2 if c in hi else 1.4)
pos={}
order=dc.iloc[-1].sort_values()
last=-9
for c,v in order.items():
    y=max(v,last+1.6); pos[c]=y; last=y
for c in dc.columns:
    ax.text(2023.8,pos[c],names[c],fontsize=7.5,color=INK if c in hi else INK2,va="center")
ax.axhspan(5,6,color=AQUA,alpha=0.18,lw=0); ax.text(1990.5,6.3,"threshold for ectomycorrhizal change in van der Linde et al. (2018): 5 to 6 kg",fontsize=7,color=INK2)
ax.set_xlim(1990,2028); ax.set_ylim(0,46); ax.set_ylabel("kg nitrogen per hectare per year")
ax.set_xticks([1990,2000,2010,2020])
ax.set_title("Nitrogen falling from the air, 1990 to 2023")
source(ax,"Total (oxidised + reduced, wet + dry) deposition, grid average, area-weighted over all administrative units. Data: EMEP MSC-W rv5.6 trend runs.")
save(fig,"ndep_countries")

r=json.load(open("out/core_nitrogen_results.json")); L=json.load(open("out/legacy_results.json"))
b=r["share_by_dep_bin"]; mids=[1.5,4,6.25,8.75,11.25,13.75,17.5,22.5,30]
sh=[v["share"] for v in b.values()]; n=[v["records"] for v in b.values()]
fig,ax=plt.subplots(figsize=(6.6,3.6))
a_,b_=L["dose_curve"]["a"],L["dose_curve"]["b"]; x=np.linspace(0.8,40,300)
ax.plot(x,1/(1+np.exp(-(a_+b_*np.log(x)))),color=BLUE,lw=2)
ax.scatter(mids,sh,s=[max(12,v/2500) for v in n],color=BLUE,edgecolor="white",linewidth=1.2,zorder=3)
ax.axvspan(5,6,color=AQUA,alpha=0.18,lw=0)
for xx,s in zip(mids,sh): ax.text(xx+0.6,s+0.025,f"{s:.0%}",fontsize=7.2,color=INK2)
ax.set_xlim(0,38); ax.set_ylim(0,1.02); ax.set_xlabel("nitrogen deposition, five-year mean (kg per hectare per year)")
ax.set_ylabel("nitrogen-sensitive share")
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
ax.set_title("Where nitrogen falls, the sensitive fungi give way")
source(ax,"Share of nitrogen-sensitive genera among indicator ectomycorrhizal records, unit-years 1990-2025, nine national recording schemes. Dot area: records.")
save(fig,"dose_curve")

T=pd.read_csv("out/legacy_trajectories.csv")
fig,ax=plt.subplots(figsize=(6.6,3.6))
for c,col in [("under 10 kg in 1990",AQUA),("10 to 20 kg",BLUE),("over 20 kg",ORANGE)]:
    t=T[(T.cls==c)&(T.ind>=100)].set_index("year").share.rolling(3,center=True,min_periods=2).mean()
    ax.plot(t.index,t.values,color=col,lw=2); ax.text(t.index[-1]+0.5,t.values[-1],c,fontsize=7.5,color=INK2,va="center")
ax.set_xlim(1990,2034); ax.set_ylim(0,1); ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
ax.set_ylabel("nitrogen-sensitive share"); ax.set_xticks([1990,2000,2010,2020])
ax.set_title("The sensitive fungi are coming back, from every starting point")
source(ax,"Units grouped by 1990 deposition. Britain, Denmark, Sweden, Belgium and Norway, long-running national schemes only. Three-year running mean.")
save(fig,"recovery_trajectories")

L=json.load(open("out/legacy_results.json"))
fig,ax=plt.subplots(figsize=(6.6,3.6))
x=np.linspace(1,30,300); ax.plot(x,1/(1+np.exp(-(a_+b_*np.log(x)))),color="#b9b8b3",lw=2,label="curve fitted across all places and years")
for c,col in [("under 10 kg in 1990",AQUA),("10 to 20 kg",BLUE),("over 20 kg",ORANGE)]:
    v=L["curve_vs_observed"][c]; (d0,d1),(s0,s1)=v["dep"],v["obs"]
    ax.annotate("",xy=(d1,s1),xytext=(d0,s0),arrowprops=dict(arrowstyle="-|>",color=col,lw=2.2,mutation_scale=14))
    ax.scatter([d0],[s0],color=col,s=40,zorder=3,edgecolor="white")
    dx,dy={"10 to 20 kg":(0.6,0.05)}.get(c,(0.5,-0.06))
    ax.text(d0+dx,s0+dy,f"{c}\n{s0:.0%} to {s1:.0%}",fontsize=7.4,color=INK2)
ax.set_xlim(0,30); ax.set_ylim(0,1.02); ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
ax.set_xlabel("nitrogen deposition (kg per hectare per year)"); ax.set_ylabel("nitrogen-sensitive share")
ax.legend(loc="upper right")
ax.set_title("Recovery ran back along the same curve")
source(ax,"Arrows: first five years (1990-94) to last five years (2021-25), long-running schemes in Britain, Denmark, Sweden, Belgium and Norway.")
save(fig,"walk_down_curve")

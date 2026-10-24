import pandas as pd, numpy as np, json
T=pd.read_parquet("data/sf_species_years.parquet")
T["fy"]=T.first_year.astype(int)
a=T.groupby("fy").size()
yrs=range(1750,2024)
a=a.reindex(yrs,fill_value=0)
dec=T.groupby((T.fy//10)*10).size()
print(dec.loc[1900:].to_string())
r=a.rolling(5).mean()
print("5y mean per yr at 1950,1980,2000,2010,2015,2020,2023:",[round(r[y]) for y in [1950,1980,2000,2010,2015,2020,2023]])
print("2014-2023 mean",a.loc[2014:2023].mean(),"2019-2023",a.loc[2019:2023].mean(), "max year",a.idxmax(),a.max())
print("share named after 2000",(T.fy>=2000).mean(),"after 1990",(T.fy>=1990).mean())
# by phylum recent
rec=T[T.fy>=2014]
print(rec.phylum.value_counts().head(8))
print(rec["class"].value_counts().head(12))
print(T.phylum.value_counts().head(8))
# time to name scenarios
cur=a.loc[2014:2023].mean()
named=len(T)
out={"named":named,"rate_2014_2023":cur,"rate_2019_2023":a.loc[2019:2023].mean()}
for est in [1.5e6,2.2e6,3.8e6,6.3e6]:
    out[f"years_to_name_{int(est)}"]=(est-named)/cur
print(out)
json.dump(out,open("out/naming_summary.json","w"),indent=1)
a.to_csv("out/descriptions_per_year.csv")

import pandas as pd, numpy as np
D="data/sf/"
N=pd.read_csv(D+"Name.tsv",sep="\t",dtype=str,quoting=3,on_bad_lines="skip")
R=pd.read_csv(D+"Reference.tsv",sep="\t",dtype=str,quoting=3,on_bad_lines="skip")
T=pd.read_csv(D+"Taxon.tsv",sep="\t",dtype=str,quoting=3,on_bad_lines="skip")
S=pd.read_csv(D+"Synonym.tsv",sep="\t",dtype=str,quoting=3,on_bad_lines="skip")
print(len(N),len(R),len(T),len(S))
print(N["rank"].value_counts().head(8))
R["year"]=pd.to_numeric(R["year"],errors="coerce")
N=N.merge(R[["ID","year"]].rename(columns={"ID":"referenceID"}),on="referenceID",how="left")
ny=N.set_index("ID")["year"]; nrank=N.set_index("ID")["rank"]
T=T[nrank.reindex(T["nameID"]).values=="sp."].copy()
T["acc_year"]=ny.reindex(T["nameID"]).values
S["year"]=ny.reindex(S["nameID"]).values
syn_min=S.groupby("taxonID")["year"].min()
T["syn_min"]=syn_min.reindex(T["ID"]).values
T["first_year"]=T[["acc_year","syn_min"]].min(axis=1)
print("accepted species",len(T), "with year",T.first_year.notna().sum())
print(T.first_year.describe())
T[["ID","phylum","class","order","family","genus","species","acc_year","first_year"]].to_parquet("data/sf_species_years.parquet")

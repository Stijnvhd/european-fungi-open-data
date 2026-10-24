import requests, json
ECM = ["Amanita muscaria","Amanita rubescens","Amanita citrina","Amanita phalloides","Boletus edulis","Imleria badia","Xerocomellus chrysenteron","Leccinum scabrum","Suillus luteus","Suillus grevillei","Paxillus involutus","Laccaria laccata","Laccaria amethystina","Russula ochroleuca","Russula nigricans","Russula cyanoxantha","Lactarius quietus","Lactarius turpis","Lactarius tabidus","Cantharellus cibarius","Hydnum repandum","Scleroderma citrinum","Tricholoma sulphureum","Cortinarius semisanguineus","Hebeloma crustuliniforme","Inocybe geophylla"]
SAP = ["Hypholoma fasciculare","Mycena galericulata","Mycena pura","Clitocybe nebularis","Lepista nuda","Gymnopus dryophilus","Rhodocollybia butyracea","Armillaria mellea","Coprinus comatus","Coprinellus micaceus","Pleurotus ostreatus","Kuehneromyces mutabilis","Phallus impudicus","Lycoperdon perlatum","Macrolepiota procera","Marasmius oreades","Agaricus campestris","Stropharia aeruginosa","Flammulina velutipes","Calocera viscosa","Mycena inclinata","Auricularia auricula-judae","Clitocybe odora","Gymnopus peronatus"]
GEN_SENS = ["Cortinarius","Tricholoma","Hydnum","Cantharellus","Hydnellum","Sarcodon","Craterellus","Phellodon"]
GEN_TOL = ["Paxillus","Laccaria","Scleroderma","Thelephora","Hebeloma"]
SP_TOL=["Lactarius quietus","Lactarius tabidus"]
out={}
def m(name,rank=None):
    p={"name":name,"kingdom":"Fungi","strict":"true"}
    if rank: p["rank"]=rank
    r=requests.get("https://api.gbif.org/v1/species/match",params=p).json()
    return r.get("usageKey"), r.get("scientificName"), r.get("matchType"), r.get("status")
for grp,lst,rank in [("ECM",ECM,"SPECIES"),("SAP",SAP,"SPECIES"),("GEN_SENS",GEN_SENS,"GENUS"),("GEN_TOL",GEN_TOL,"GENUS"),("SP_TOL",SP_TOL,"SPECIES")]:
    for n in lst:
        k,s,mt,st=m(n,rank); out.setdefault(grp,[]).append(dict(name=n,key=k,sci=s,match=mt,status=st)); print(grp,n,k,mt,st)
json.dump(out,open("data/keys.json","w"),indent=1)

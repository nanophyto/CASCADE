import pandas as pd
import numpy as np
from yaml import  load, Loader
from pandas.api.types import is_integer_dtype
from datetime import datetime


with open('/home/phyto/CoccoData/phases.yml', 'r') as f:
    phases = load(f, Loader=Loader)


with open('/home/phyto/CoccoData/family.yml', 'r') as f:
    families = load(f, Loader=Loader)


d = pd.DataFrame.from_dict(phases, orient='index')
d = d.rename_axis("species").reset_index()
d['genera'] = d['species'].str.split(" ").str[0]

inverse = {}
for k,v in families.items():
    for x in v:
        inverse.setdefault(x, []).append(k)

df = pd.DataFrame.from_dict(inverse, orient='index')
df = df.rename_axis("genera").reset_index()
df = df.rename(columns={0: "family"})

d = pd.merge(d, df, on='genera', how="outer")

d.to_csv("/home/phyto/CoccoData/library.csv", index=False)



print("fin")
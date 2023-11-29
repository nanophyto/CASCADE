
import pandas as pd
import numpy as np
import sys
import geopandas as gpd
from shapely.geometry import Point # Point class

from yaml import safe_load, load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

import yaml



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

classification_dict = d.set_index('species').to_dict('index')



def import_villiot2021():    
    #read obrien size data
    villiot_size = pd.read_csv("/home/phyto/CoccoData/viliot2021_cell_diameters.csv")

    #resample ehux
    RCC1731 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1731"]['mean'], villiot_size[villiot_size["strain"]=="RCC1731"]['std'], 10000)
    RCC1128 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1228"]['mean'], villiot_size[villiot_size["strain"]=="RCC1228"]['std'], 10000)
    
    ehux = {'species': ["Emiliania huxleyi"],
                'mean': [np.mean([RCC1128, RCC1731])],
                'std': [np.std([RCC1128, RCC1731])],
                'strain': ['NA']}
    ehux = pd.DataFrame(ehux)
    villiot_size = villiot_size[villiot_size["species"]!="Emiliania huxleyi"]
    
    #resample leptoporus
    RCC1130 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1130"]['mean'], villiot_size[villiot_size["strain"]=="RCC1130"]['std'], 10000)
    RCC1135 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1135"]['mean'], villiot_size[villiot_size["strain"]=="RCC1135"]['std'], 10000)
    
    leptoporus = {'species': ["Coccolithus leptoporus"],
                'mean': [np.mean([RCC1130, RCC1135])],
                'std': [np.std([RCC1130, RCC1135])],
                'strain':['NA']}
    leptoporus = pd.DataFrame(leptoporus)

    villiot_size = villiot_size[villiot_size["species"]!="Calcidiscus leptoporus"]

    villiot_size = pd.concat([villiot_size, ehux, leptoporus])

    villiot_size = villiot_size.rename(columns={'mean': "size_mean_villiot2021",
                                          'std': "size_std_villiot2021"})

    villiot2021 = villiot_size[['species', 'size_mean_villiot2021', 'size_std_villiot2021']]

    villiot_dict = villiot2021.set_index('species').to_dict('index')

    return(villiot_dict)

villiot2021 = import_villiot2021()




def import_obrien2013():
    
    #read obrien size data
    obrien_size = pd.read_csv("/home/phyto/CoccoData/obrien_cell_diameters.csv")
    #apply synonyms
    obrien_size['species'] = obrien_size['species'].str.strip()

    with open('/home/phyto/CoccoData/groupings.yml', 'r') as f:
        groupings = load(f, Loader=Loader)

    species = obrien_size['species']
    dict = {species:k
        for k, v in groupings.items()
        for species in v['alt']}

    obrien_size = obrien_size.replace(dict)

    df = obrien_size.rename(columns={'diameter': "size_mean_obrien2013"})

    #take mean for duplicate entries
    df = df.groupby('species').mean().reset_index()    

    ##########
    # NEED TO ESTIMATE VOLUME FROM DIAMETER
    ##########

    obrien_dict = df.set_index('species').to_dict('index')

    #merge size data
    return(obrien_dict)

obrien2013 = import_obrien2013()




def merge(a: dict, b: dict, path=[]):
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] != b[key]:
                raise Exception('Conflict at ' + '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a


dict1 = merge(villiot2021, obrien2013)


dict2 = merge(dict1, classification_dict)


with open('/home/phyto/CoccoData/test.yml', 'w') as outfile:
    yaml.dump(dict2, outfile, default_flow_style=False)


print("fin")
import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader

def rename_synonyms(d, index='species', remove_duplicate=True):
    d['species'] = d['species'].str.strip()

    with open('/home/phyto/CoccoData/data/classification/synonyms.yml', 'r') as f:
        groupings = load(f, Loader=Loader)

    synonym_dict = {species:k
        for k, v in groupings.items()
        for species in v}

    d_new = d.replace(synonym_dict)

    species_original = d['species']
    species_renamed = d_new['species']

    changed_species = pd.DataFrame(set(species_original).difference(species_renamed))

    renamed_species_renamed =  changed_species.replace(synonym_dict)

    d = pd.concat([changed_species, renamed_species_renamed], axis=1)

    print(d)

    d.to_csv("/home/phyto/CoccoData/data/unprocessed/sheward2024_30nov_renamed.csv")

    return(d)


d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sheward2024_30nov.csv")
d['species'] = d['genus'] + " " + d['Species']
rename_synonyms(d)

print("fin")
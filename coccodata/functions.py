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

    d = d.replace(synonym_dict)

    #take mean for duplicate entries
    if remove_duplicate == True: 
        d = d.groupby(index).mean().reset_index()    



    #check if final species are in species library
    species_library = {v: k for k, v in synonym_dict.items()}

    species_observed = d['species']

    if (set(species_observed).issubset(species_library)):
        print("all species are defined")
    else:
        raise ValueError("undefined species:" + str(set(species_observed).difference(species_library)))

    return(d)
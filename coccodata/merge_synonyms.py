


import pandas as pd
import numpy as np
from yaml import load, Loader
from yaml import dump, Dumper

with open('/home/phyto/CoccoData/data/unprocessed/synonyms/synonyms.yml', 'r') as f:
            synonyms1 = load(f, Loader=Loader)

with open('/home/phyto/CoccoData/data/unprocessed/synonyms/synonyms2.yml', 'r') as f:
            synonyms2 = load(f, Loader=Loader)

s = synonyms1 | synonyms2

path = '/home/phyto/CoccoData/data/classification/synonyms.yml'
with open(path, 'w') as outfile:
    dump(s, outfile, default_flow_style=False)

print("fin")
# def merge_synonyms():

#     with open('/home/phyto/CoccoData/data/unprocessed/synonyms/synonyms.yml', 'r') as f:
#                 synonyms1 = load(f, Loader=Loader)

#     with open('/home/phyto/CoccoData/data/unprocessed/synonyms/synonyms2.yml', 'r') as f:
#                 synonyms2 = load(f, Loader=Loader)


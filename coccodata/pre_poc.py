import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 
import glob, os
sys.path.insert(0, '/home/phyto/CoccoData/coccodata/')
from functions import rename_synonyms



def import_data(path):

    all_files = glob.glob(os.path.join(path, "*.csv"))

    d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
    d = d.fillna(0)

    #rename synonyms and typos:
    d = rename_synonyms(d, remove_duplicate=False, check_synonyms = False)

    return(d)




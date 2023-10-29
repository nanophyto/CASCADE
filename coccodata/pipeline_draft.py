import pandas as pandas
import numpy as np 
from collections import namedtuple




# pipeline for size estinates:
lab = [True, False]
number_of_lab_studies = [1, 3]
lab_std = [True, False]
species = ["spp_a", "spp_b"]

hol = [True, False]
alternative_phase = [None, "spp_a"]

amt = [True, True]
amt_std = [True, True]

Nannotax = [True, True]



for i in range(len(species)):

    if lab[i]:
        if number_of_lab_studies==1:
            if lab_std ==True:
                print("resample from study")
            else:
                print("use mean")
        else:
            print("nested resampling with mean or std")
    
    elif amt[i]:
        if amt_std ==True:
            print("resample from study")
        else:
            print("use mean")

    elif (hol[i] & alternative_phase[i]!=None):
        print("use alternative phase")

    elif Nannotax[i]:
        print("use NannoTax data")

    else:
        raise ValueError("size not defined for any method")



"""

Pipeline for POC:

Loop through species and estimate POC based on size and allometric scaling.

Both have uncertainties and are resampled

"""

sizes_dict = {"Coccolithus leptoporus": {
            'mean': [10],
            'std': [1]},
            "Emiliania huxleyi": {
            'mean': [10],
            'std': [1]}
}

sample_n = 10000


allometry: {
    'a':{
        'mean': [0.8],
        'std': [0.1]},
    'b':{
        'mean': [0.3],
        'std': [0.04]}
}
    
POC_dict = dict()

for i in range(len(species)):

    sp = species[i]
    sizes_list = np.random.normal(sizes_dict[sp]['mean'], sizes_dict[sp]['std'], sample_n)
    a_list = np.random.normal(allometry['a']['mean'], allometry['a']['std'], sample_n)
    b_list = np.random.normal(allometry['b']['mean'], allometry['b']['std'], sample_n)

    POC_list = a_list*sizes_list**b_list

    POC_mean = np.mean(POC_list)
    POC_std = np.std(POC_list)

    print("append new values to POC dictionary")


"""

Pipeline for PIC:


Estimate PIC for species where availible

Estimate PIC:POC for:
    - HOL 
    - Genera
    - Families

    
Estimate PIC for each species by multiplying species POC
with group PIC:POC.

With the following order:
    1. HOL
    2. Genera
    3. Family

"""

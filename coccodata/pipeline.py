import pandas as pandas
import numpy as np 
from collections import namedtuple
from yaml import load, Loader
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/home/phyto/CoccoData/coccodata/')
from library import library
n = 1000

m = library('/home/phyto/CoccoData/data/classification/phases.yml',
            '/home/phyto/CoccoData/data/classification/family.yml',
            "/home/phyto/CoccoData/data/sizes/")

#return the named tuple so we can use it:
ntpl = m.return_ntpl()
#also return species list:
species_list = m.return_species_list()
name = species_list[1]
#name = "Helicosphaera pavimentum HOL"

def check_HOL(ntpl, name):
    """
    returns name of HET phase 
    if size data of HOL is undefined
    and if alternative phase is defined
    """

    species_library =  [t for t in ntpl  if t.species == name]
    if all(d.mean is None for d in species_library): 
        if (any(d.phase == "HOL" for d in species_library)):
            if (any(d.alternate_phase is not None for d in species_library)):
                return(species_library[0].alternate_phase)
    else:
        return(None)

def resample_sd(size_library):
    #estimate mean SD across studies by bootstrapping
    sd_library = np.array([ x.sd for x in size_library])
    sd_library = sd_library[sd_library != np.array(None)]
    sd_library = sd_library[sd_library > 0]
    sd_estimate = np.random.choice(sd_library, n)
    return(sd_estimate)
    
def resample_size(ntpl, spp_name):

    if check_HOL(ntpl, spp_name) !=None:
        HET_name = check_HOL(ntpl, spp_name)
        size_library =  [t for t in ntpl  if t.species == HET_name]
        print(HET_name)
    else:
        size_library =  [t for t in ntpl  if t.species == spp_name]

    #simulate SDs for studies which do not have a SD:
    sd_estimate = resample_sd(size_library)

    #simulate size distributions for each study which has data: 
    estimate = []

    for i in range(len(size_library)):

        if (size_library[i].mean is not None) and ((size_library[i].sd is None) or (size_library[i].sd ==0)):
            size_estimate = []
            for j in range(n):
                size_estimate_n = np.random.normal(size_library[i].mean, sd_estimate[j], 1)
                size_estimate.extend(size_estimate_n)

        elif (size_library[i].mean is None) and (size_library[i].sd is None):
            size_estimate = []
        else:
            size_estimate = np.random.normal(size_library[i].mean, size_library[i].sd, n)
            print(size_library[i].id)
        size_estimate = [x for x in size_estimate if x >= 0 ]
        estimate.extend(size_estimate)

    return(estimate)

estimate = resample_size(ntpl, name)
sns.histplot(x=estimate)
plt.show()


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")

import pandas as pd
d = pd.DataFrame({'spp1':estimate, 'spp2': estimate})
d = pd.melt(d, id_vars=['sp'], value_vars=['B'])


N_plots_per_page = 9
for cols in grouper(d['species'].unique(), N_plots_per_page):
    g = sns.catplot(data=d, x='solution', y='score', col='subject', col_wrap=3, kind='point', col_order=cols)
    pdf.savefig(g.fig)
pdf.close()

spp_list = []

for i in range(len(species_list)):
    name = species_list[i]
    species_library =  [t for t in library  if t.species == name]

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

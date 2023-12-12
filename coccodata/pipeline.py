import pandas as pandas
import numpy as np 
from collections import namedtuple
from yaml import load, Loader
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.insert(0, '/home/phyto/CoccoData/coccodata/')
from library import library
from regression import regression_simulation
from plotting import multipage_plots
n = 1000

m = library('/home/phyto/CoccoData/data/classification/phases.yml',
            '/home/phyto/CoccoData/data/classification/family.yml',
            "/home/phyto/CoccoData/data/sizes/")

#return the named tuple so we can use it:
ntpl = m.return_ntpl()
#also return species list:
species_list = m.return_species_list()
name = species_list[0]
#name = "Helicosphaera pavimentum HOL"
#name = "Emiliania huxleyi"


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
    if sd_library.any():
        sd_estimate = np.random.choice(sd_library, n)
    else:
        sd_estimate = [0]*n
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

def resample_POC(size_simulation):
    poc_simulation = []
    
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
    Y_train = d['pg poc']
    X_train = d['volume']

    m = regression_simulation(X_train, size_simulation, Y_train)

    poc_simulation = m.simulate_data()

    return(poc_simulation)


def bootstrap_measurements(species_name):

    size_simulation = resample_size(ntpl, species_name)
    poc_simulation = resample_POC(size_simulation)
    diameter_simulation = (6*np.asarray(size_simulation)/np.pi)**(1/3)


    d_vol = pd.DataFrame({'species': species_name, 'value': size_simulation, 'variable':'volume'})
    d_poc = pd.DataFrame({'species': species_name, 'value': poc_simulation, 'variable':'pg poc'})
    d_dia = pd.DataFrame({'species': species_name, 'value': diameter_simulation, 'variable':'diameter'})

    d = pd.concat([d_dia, d_vol, d_poc])

    return(d)

estimates = []

for i in range(len(species_list)): #
    print(species_list[i])
    estimates.append(bootstrap_measurements(species_list[i]))

estimates = pd.concat(estimates)

#estimates.to_csv("/home/phyto/CoccoData/test.csv")

#estimates = pd.read_csv("/home/phyto/CoccoData/test.csv")
multipage_plots(d=estimates, species_list=species_list, n_page=8, out_path='/home/phyto/multipage_pdf.pdf')

print("fin")
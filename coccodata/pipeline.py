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

n = 10000

d = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")
species_list = d['species']

m = library('/home/phyto/CoccoData/data/classification/phases.yml',
            '/home/phyto/CoccoData/data/classification/family.yml',
            "/home/phyto/CoccoData/data/sizes/",
            "/home/phyto/CoccoData/data/pic/",
            species_list)

#create a list of HOL species for future use:
HOL_list = m.return_HOL_list()

#return the named tuple so we can use it:
ntpl = m.return_ntpl()
ntpl_size =  [t for t in ntpl  if t.measurement == "volume"]
ntpl_pic =  [t for t in ntpl  if t.measurement == "pic"]
ntpl_poc =  [t for t in ntpl  if t.measurement == "poc"]

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
    print("estimating cross study sd for spp")
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

    print("resampling size for spp: " + spp_name)

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
        #    print(size_library[i].id)
        size_estimate = [x for x in size_estimate if x >= 0 ]
        estimate.extend(size_estimate)

    return(estimate)


def c_regression(measurement, size_simulation):
    #function to estimate carbon based on size distribution

    if measurement == "pic":
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
        Y_train = d['pg poc']
        X_train = d['volume']
    elif measurement == "poc":
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
        Y_train = d['pg poc']
        X_train = d['volume']
    else:
        raise ValueError("measurement undefined, should be 'pic' or 'poc'")

    m = regression_simulation(X_train, size_simulation, Y_train)
    # m.return_performance()
    c_simulation = m.simulate_data()

    return(c_simulation)


    
def resample_carbon_species(ntpl, spp_name, n, size_simulation, measurement):

    library =  [t for t in ntpl  if t.species == spp_name]

    print("estimating " + measurement + " for spp " + spp_name)

    #simulate size distributions for each study which has data: 
    estimate = []

    #drop empty studies from library
    library =  [t for t in library  if t.mean != None]

    if len(library)>0:
        try:
            #simulate SDs of all species to be used for studies which do not have a SD:
            cross_study_sd_estimate = resample_sd(library)
        #if direct measurements exist use them:
            print("estimating "+ measurement + " based on resampling")
            for i in range(len(library)):
                #if mean is known, bootstrap:
                if (library[i].sd is None) or (library[i].sd ==0):
                    #if sd is unknown, use cross_study resampled sd
                    estimate.extend(np.random.normal(library[i].mean, cross_study_sd_estimate[i], n))
                else:
                    #if sd is known, use study sd
                    estimate.extend( np.random.normal(library[i].mean, cross_study_sd_estimate[i], n))
        except:
            #if resample of sd fails because there are no SDs
            #estimate carbon based on size:
            print("estimating "+ measurement + " based on GLM")
            estimate.extend(c_regression(measurement, size_simulation))                
    else:
    #otherwise estimate carbon based on size:
        print("estimating "+ measurement + " based on GLM")
        estimate.extend(c_regression(measurement, size_simulation))    

    return(estimate)



def resample_measurements(species_name):

    size_simulation = resample_size(ntpl_size, species_name) 
    diameter_simulation = (6*np.asarray(size_simulation)/np.pi)**(1/3)

    pic_simulation = resample_carbon_species(ntpl_pic, species_name, n, size_simulation, "pic")
    poc_simulation = resample_carbon_species(ntpl_poc, species_name, n, size_simulation, "poc")

    d_vol = pd.DataFrame({'species': species_name, 'value': size_simulation, 'variable':'volume'})
    d_poc = pd.DataFrame({'species': species_name, 'value': poc_simulation, 'variable':'pg poc'})
    d_dia = pd.DataFrame({'species': species_name, 'value': diameter_simulation, 'variable':'diameter'})
    d_pic = pd.DataFrame({'species': species_name, 'value': pic_simulation, 'variable':'pg pic'})

    d = pd.concat([d_dia, d_vol, d_poc, d_pic])
    #d = pd.concat([d_dia, d_vol, d_poc])

    return(d)


#THIS IS CURRENTLY SUPER SLOW!
estimates = []

for i in range(len(species_list)): #
     print(species_list[i])
     estimates.append(resample_measurements(species_list[i]))

# for i in range(8): #
#    print(species_list[i])
#    estimates.append(resample_measurements(species_list[i]))

estimates = pd.concat(estimates)
estimates.to_csv("/home/phyto/CoccoData/cellular_dataset.csv", index=False)

d = pd.read_csv("/home/phyto/CoccoData/cellular_dataset.csv")

# species_list = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")
# species_list = list(species_list['species'])

# estimate = estimate_group_pic_poc(d, species_list, 
#             '/home/phyto/CoccoData/data/classification/phases.yml',
#             '/home/phyto/CoccoData/data/classification/family.yml')

# estimate.append_PIC_list(group="HOL")
# d = estimate.return_d()      
# estimate.append_PIC_list(group="genera")
# d = estimate.return_d()
# estimate.append_PIC_list(group="family")
# d = estimate.return_d()

# d.to_csv("/home/phyto/CoccoData/final_pic.csv")

multipage_plots(d=d, species_list=species_list[0:7], n_page=8, out_path='/home/phyto/multipage_pdf.pdf')

#estimates.to_csv("/home/phyto/CoccoData/test.csv")

#estimates = pd.read_csv("/home/phyto/CoccoData/test.csv")
#multipage_plots(d=d, species_list=species_list, n_page=8, out_path='/home/phyto/library.pdf')

print("fin")
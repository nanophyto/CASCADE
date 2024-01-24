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
from functions import bootstrap

class pipeline:

    def __init__(self, root= '/home/phyto/CoccoData/data/', n = 1000):

        d = pd.read_csv(root + "species_list.csv")
        self.species_list = d['species']

        m = library(root + 'classification/phases.yml',
                    root + 'classification/family.yml',
                    root + "sizes/",
                    root + "pic/",
                    root + "poc/",           
                    self.species_list)

        #create a list of HOL species for future use:
        self.HOL_list = m.return_HOL_list()

        #return the named tuple so we can use it:
        self.ntpl = m.return_ntpl()
        self.ntpl_size =  [t for t in self.ntpl  if t.measurement == "volume"]
        self.ntpl_pic =  [t for t in self.ntpl  if t.measurement == "pic"]
        self.ntpl_poc =  [t for t in self.ntpl  if t.measurement == "poc"]

        self.n = n 
    def return_species_list(self):
        return(self.species_list)

    def check_HOL(self, name):
        """
        returns name of HET phase 
        if size data of HOL is undefined
        and if alternative phase is defined
        """
        species_library =  [t for t in self.ntpl  if t.species == name]
        if all(d.mean is None for d in species_library): 
            if (any(d.phase == "HOL" for d in species_library)):
                if (any(d.alternate_phase is not None for d in species_library)):
                    return(species_library[0].alternate_phase)
        else:
            return(None)

    def resample_cv(self, size_library):
        print("estimating cross study cv for spp")
        #estimate mean SD across studies by bootstrapping
        sd_library = np.array([ x.sd for x in size_library if (x.sd is not None) and (x.mean is not None)])
        mean_library = np.array([x.mean for x in size_library if (x.sd is not None) and (x.mean is not None) ])
        #drop entries where sd values are None
        sd_library = sd_library[sd_library != np.array(None)]
        mean_library = mean_library[sd_library > 0]
        sd_library = sd_library[sd_library > 0]
        #convert sd to cv by dividing by mean:
        cv_library = sd_library/mean_library

        #bootstrap cv estimates
        if cv_library.any():
            cv_estimate = bootstrap(cv_library, self.n)
        else:
            cv_estimate = [0]*self.n

        return(cv_estimate)
        
    def resample_size(self, spp_name):

        print("resampling size for spp: " + spp_name)

        if self.check_HOL(spp_name) !=None:
            HET_name = self.check_HOL(self.ntpl_size, spp_name)
            size_library =  [t for t in self.ntpl  if t.species == HET_name]
            print(HET_name)
        else:
            size_library =  [t for t in self.ntpl  if t.species == spp_name]

        #simulate SDs for studies which do not have a SD:
        all_studies_cv = self.resample_cv(size_library)
        all_species_cv = self.resample_cv([t for t in self.ntpl])

        #simulate size distributions for each study which has data: 
        estimate = []

        for i in range(len(size_library)):
            #for every study:
            if (size_library[i].mean is not None) and ((size_library[i].sd is None) or (size_library[i].sd ==0)):
                size_estimate = []
                if np.mean(all_studies_cv)  > 0:
                    for j in range(self.n):
                        random_sd =  all_studies_cv[j]*size_library[i].mean
                        size_estimate_n = np.random.normal(size_library[i].mean, random_sd, 1)
                        size_estimate.extend(size_estimate_n)
                else:
                    print("group sd unknown, using cross species sd instead (size)")
                    for j in range(self.n):
                        #if 0 or NA:
                        random_sd =  all_species_cv[j]*size_library[i].mean
                        size_estimate_n = np.random.normal(size_library[i].mean, random_sd, 1)
                        size_estimate.extend(size_estimate_n)

            elif (size_library[i].mean is None) and (size_library[i].sd is None):
                size_estimate = []
            else:
                size_estimate = np.random.normal(size_library[i].mean, size_library[i].sd, self.n)
            #    print(size_library[i].id)
            size_estimate = [x for x in size_estimate if x > 0 ]
            estimate.extend(size_estimate)

        return(estimate)


    def c_regression(self, name, measurement, size_simulation):
        #function to estimate carbon based on size distribution

        if measurement == "pic":
            d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/rosie_size_pic.csv")
            d = d.dropna()
            d = d[d['PIC pg C'] >0]
            d = d[d['Volume'] >0]
            Y_train = d['PIC pg C']
            X_train = d['Volume']

            species_library =  [t for t in self.ntpl_pic  if t.species == name]


            #TO ADD ONE HOT ENCODED PHASES!!!

        elif measurement == "poc":
            d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
            Y_train = d['pg poc']
            X_train = d['volume']
        else:
            raise ValueError("measurement undefined, should be 'pic' or 'poc'")

        r = regression_simulation(X_train, size_simulation, Y_train)
        # m.return_performance()
        c_simulation = r.simulate_data()
        return(c_simulation)


        
    def resample_carbon_species(self, spp_name, size_simulation, measurement):

        library =  [t for t in self.ntpl  if (t.species == spp_name and
                                            t.measurement == measurement)]

        print("estimating " + measurement + " for spp " + spp_name)

        #simulate size distributions for each study which has data: 
        estimate = []

        #drop empty studies from library
        library =  [t for t in library  if t.mean != None]

        if len(library)>0:
            try:
                #simulate SDs of all species to be used for studies which do not have a SD:
                all_studies_cv = self.resample_cv(library)
                all_species_cv = self.resample_cv([t for t in self.ntpl])

            #if direct measurements exist use them:
                print("estimating "+ measurement + " based on resampling")
                for i in range(len(library)):
                    #if mean is known:
                    if (library[i].sd is None) or (library[i].sd ==0):
                        #if sd is not known
                        if np.mean(all_studies_cv)  > 0:
                            #if species sd is known
                            for j in range(self.n):
                                random_sd =  all_studies_cv[j]*library[i].mean
                                estimate.extend(np.random.normal(library[i].mean, random_sd, 1))
                        else:
                            #if species sd is unknown
                            print("group sd unknown, using cross species sd instead")
                            for j in range(self.n):
                                #if 0 or NA:
                                random_sd =  all_species_cv[j]*library[i].mean
                                estimate.extend(np.random.normal(library[i].mean, random_sd, 1))
                    else:
                        #if sd is known, use study sd
                        estimate.extend( np.random.normal(library[i].mean, library[i].sd, self.n))
            except:
                #if resample of sd fails because there are no SDs
                #estimate carbon based on size:
                print("estimating "+ measurement + " based on GLM")
                estimate.extend(self.c_regression(spp_name, measurement, size_simulation))                
        else:
        #otherwise estimate carbon based on size:
            print("estimating "+ measurement + " based on GLM")
            estimate.extend(self.c_regression(spp_name, measurement, size_simulation))    

        return(estimate)



    def resample_measurements(self, species_name):

        size_simulation = self.resample_size(species_name) 
        diameter_simulation = (6*np.asarray(size_simulation)/np.pi)**(1/3)

        pic_simulation = self.resample_carbon_species(species_name, size_simulation, "pic")
        poc_simulation = self.resample_carbon_species(species_name, size_simulation, "poc")

        d_vol = pd.DataFrame({'species': species_name, 'value': size_simulation, 'variable':'volume'})
        d_poc = pd.DataFrame({'species': species_name, 'value': poc_simulation, 'variable':'pg poc'})
        d_dia = pd.DataFrame({'species': species_name, 'value': diameter_simulation, 'variable':'diameter'})
        d_pic = pd.DataFrame({'species': species_name, 'value': pic_simulation, 'variable':'pg pic'})

        d = pd.concat([d_dia, d_vol, d_poc, d_pic])
        #d = pd.concat([d_dia, d_vol, d_poc])

        return(d)

#test = resample_measurements("Gephyrocapsa ericsonii")

m = pipeline()
species_list = m.return_species_list()
#THIS IS CURRENTLY SUPER SLOW!
estimates = []

for i in range(len(species_list)): #
     print(species_list[i])
     estimates.append(m.resample_measurements(species_list[i]))
     print("finished estimating species #" + str(i+1) + " out of " + str(len(species_list)) + " species")

# for i in range(6): #
#      print(species_list[i])
#      estimates.append(m.resample_measurements(species_list[i]))

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

#multipage_plots(d=d, species_list=species_list[0:7], n_page=8, out_path='/home/phyto/multipage_pdf.pdf')

m = multipage_plots(d=d, species_list=species_list, n_page=6, out_path='/home/phyto/library.pdf')
m.export_pdf()

#estimates.to_csv("/home/phyto/CoccoData/test.csv")

#estimates = pd.read_csv("/home/phyto/CoccoData/test.csv")
#multipage_plots(d=d, species_list=species_list, n_page=8, out_path='/home/phyto/library.pdf')

print("fin")
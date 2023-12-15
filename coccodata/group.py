import pandas as pandas
import numpy as np 
from collections import namedtuple
from yaml import load, Loader
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.insert(0, '/home/phyto/CoccoData/coccodata/')
n = 10000

class group_library():
    '''
    A function to 

    Parameters
    ----------
        d : DataFrame
            pic and poc data

    Returns
    -------

    Example
    -------

    '''
    def __init__(self, species_list, family_yml, phases_yml,):
        self.species_list = species_list
        with open(family_yml, 'r') as f:
            self.families = load(f, Loader=Loader)


    def genera_library(self, species):
        '''
        A function to create a list of species for each genera in 
        a species list.

        ''' 
        genera_list = [x.split(' ')[0] for x in self.species_list]

        d = pd.DataFrame({'species':self.species_list, 'genera':genera_list})

        genera = species.split(" ")[0]
        values = list(d[d['genera']==genera]['species'])

        #DROP HOL!!!
        #HOL_list = ['']
        #genera_list = genera_list[~ species_list.species.isin(HOL_list)] #drop HOL from df

        return(values)
    

    def family_library(self, species):
        
        species_family = {}
        for k,v in self.families.items():
            for x in v:
                species_family.setdefault(x, []).append(k)

        genera = species.split(" ")[0]
        family = self.species_family[genera]  
        family_species = self.familiess[family[0]]  
      
        #DROP HOL!!!
        #HOL_list = ['']
        #genera_list = genera_list[~ species_list.species.isin(HOL_list)] #drop HOL from df

        return(family_species) #should return list of species in family
    
    def library(self, species, group):
        if group=="family":
            l = self.family_library(species)
        elif group == "genera":
            l = self.genera_library(species) 
        else:
            raise ValueError("group should be 'genera' or 'family'")
        return(l) # list of all species in group


'''

input:
phases.yml
family.yml
species_list = list of all species in study

output:


example:

d = pd.read_csv("/home/phyto/CoccoData/species_pic.csv")
species_list = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")

estimate = estimate_group_pic_poc(d, species_list, 
            /home/phyto/CoccoData/data/classification/phases.yml',
            '/home/phyto/CoccoData/data/classification/family.yml')
            
estimate.append_PIC_list(group="genera")
d = m.return_d()
estimate.append_PIC_list(group="family")
d = m.return_d()


'''

class estimate_group_pic_poc():
    '''
    A function to simulate PIC:POC estimates of a genera
    based on previously computed PIC and POC values.

    Parameters
    ----------
        d : DataFrame
            pic and poc data
        species_list : DataFrame
            dataframe of species list
        n: int
            Number of random samples to simulate

    Returns
    -------
    array of PIC:POC values of length n

    Example
    -------
    
    d = pd.read_csv("/home/phyto/CoccoData/test_pic.csv")
    yml_path = '/home/phyto/CoccoData/data/classification/family.yml'
    m = estimate_group_pic_poc(d, 1000, yml_path)
    estimate = m.estimate("Emiliania huxleyi")
    print(np.median(estimate))
    '''
    def __init__(self, d, species_list, family_yml, phases_yml,  n=10000):
        d['genera'] = d['species'].str.split(" ").str[0]
        self.d = d
        self.n = n
        self.species_list =  species_list
        self.grouper = group_library(phases_yml, family_yml,
            species_list)
        

    def merge_PIC_POC(self, species_name):
        d = self.d
        pic_poc = []
        #determine group:
        group_list = self.grouper.library(species_name, self.group)

        for i in range(len(group_list)):
            df = d[d['species']==group_list[i]]
            pic = list(df[df['variable']=="pg pic"]['value'])
            poc = list(df[df['variable']=="pg poc"]['value'])
            poc = [num for num in poc if num > 0] #drop negative poc values
            pic = [num for num in pic if num > 0] #drop negative poc values

            if pic and poc:
                pic =  np.random.choice(pic, self.n, replace=True)
                poc =  np.random.choice(poc, self.n, replace=True)
                pic_poc.extend(pic/poc)
            else:
                None
                
        if pic_poc:
            pic_poc =  np.random.choice(pic_poc, n, replace=True)

        return(pic_poc)


    def resample_PIC_POC(self, species_name):
        pic_poc_simulation = self.merge_PIC_POC(species_name)
        d = pd.DataFrame({'species': species_name, 'value': pic_poc_simulation, 'variable':'pic/poc'})
        return(d)


    def resample_PIC(self, species_name):
        df = self.d[self.d['species']==species_name]
        poc = list(df[df['variable']=="pg poc"]['value'])
        pic_poc = list(df[df['variable']=="pic/poc"]['value'])

        if not pic_poc:
            print("pic/poc undefined")
            print(species_name)
        if not poc:
            print("poc undefined")
            print(species_name) 
                
        if pic_poc and poc:
            pic_poc =  np.random.choice(pic_poc, self.n, replace=True)
            poc =  np.random.choice(poc, self.n, replace=True)
            try:
                pic = pic_poc*poc
            except:
                print("there is an error pic*pic_poc")
                print(species_name)
            d = pd.DataFrame({'species': species_name, 'value': pic, 'variable':'pg pic'})
            return(d)
        else:
            #print("fml")
            return(None)

    def append_PIC_POC(self, species_name):
        
        #find species for which PIC is already defined:
        pic_defined_spp = list(self.d[self.d['variable']=="pg pic"]['species'])

        if species_name in pic_defined_spp:
            None
            #print("pic of species already defined")
            #print(species_name)
        else:
            try:
                estimated_pic_poc = self.resample_PIC_POC(species_name)
                return(estimated_pic_poc)
            except:
                print("species not HOL and genera not defined")
                print("(" + species_name + ")")
                return(None)       


    def append_PIC(self, species_name):        
        #find species for which PIC is already defined:
        pic_defined_spp = list(self.d[self.d['variable']=="pg pic"]['species'])

        if species_name in pic_defined_spp:
            #print("pic already defined")
            #print(species_name)
            None
        else:
            try:
                estimated_pic = self.resample_PIC(species_name)
                return(estimated_pic)
            except:
                None
                #print("species not HOL and family not defined")
                print("(" + species_name + ")")
                return(None)


    def append_PIC_list(self, group="genera"):
        self.group = group
        list_of_species = list(self.species_list['species'])
        #append PIC:POC values based on genera:
        for i in range(len(list_of_species )): #
            new_data = self.append_PIC_POC(list_of_species[i])
            self.d = pd.concat([self.d, new_data], ignore_index=True)
        #then estimate PIC:
        for i in range(len(list_of_species)): #
            new_data = self.append_PIC(list_of_species[i])
            self.d = pd.concat([self.d, new_data], ignore_index=True)


    def return_d(self):
        return(self.d)

# d = pd.read_csv("/home/phyto/CoccoData/species_pic.csv")
# species_list = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")

# estimate = estimate_group_pic_poc(d, species_list, 
#             '/home/phyto/CoccoData/data/classification/phases.yml',
#             '/home/phyto/CoccoData/data/classification/family.yml')
            
# estimate.append_PIC_list(group="genera")
# d = estimate.return_d()
# estimate.append_PIC_list(group="family")
# d = estimate.return_d()
# print("fin")
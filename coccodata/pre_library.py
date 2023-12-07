import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 

def rename_synonyms(d, index='species', remove_duplicate=True):
    d['species'] = d['species'].str.strip()

    with open('/home/phyto/CoccoData/data/classification/synonyms.yml', 'r') as f:
        groupings = load(f, Loader=Loader)

    species = d['species']
    synonym_dict = {species:k
        for k, v in groupings.items()
        for species in v}

    d = d.replace(synonym_dict)

    #take mean for duplicate entries
    if remove_duplicate == True: 
        d = d.groupby(index).mean().reset_index()    

    return(d)

def resample_size(d, species):
    d = d[d["species"]==species]

    size = []
    for i in range(len(d)):
        size_distribution = np.random.normal(d.iloc[i]['mean'], d.iloc[i]['sd'], 10000)
        size.append(size_distribution)

    print(size)

    species_data = {'species': [species],
                'mean': [np.mean(size)],
                'sd': [np.std(size)]}
    
    species_data = pd.DataFrame(species_data)

    return(species_data)

d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/viliot2021_cell_diameters.csv")

def resample_size_d(d):
    species_list = d['species'].unique()

    new_d = []

    for i in range(len(species_list)):
        size = resample_size(d, species_list[i])
        new_d.append(size)

    d = pd.concat(new_d)

    return(d)

def pre_villiot2021a():    
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/viliot2021_cell_diameters.csv")
    d = rename_synonyms(d, ['species', 'strain'])
    d = d.rename(columns={'std': "sd"})

    d = resample_size_d(d)
    d['reference'] = "villiot2021a"
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/viliot2021a.csv", index=False)

def pre_villiot2021b():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/villiot2021_literature_morphometric_volume.csv")
    d = rename_synonyms(d)
    d = d.rename(columns={'cell volume': "mean"})
    d['mean'] = np.round(d['mean'])
    d['reference'] = "villiot2021b"
    d['method'] = 'literature morphometrics'
    d['sd'] = None
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/viliot2021b.csv", index=False)

def pre_obrien2013a():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/obrien_cell_diameters.csv")
    d = rename_synonyms(d)
    d = d.rename(columns={'diameter': "mean"})
    d['mean'] = (1/6)*math.pi*(d['mean']**3)
    d['mean'] = np.round(d['mean'])
    d['reference'] = 'obrien2013a'
    d['sd'] = None
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/obrien2013a.csv", index=False)
    return(d)

def pre_obrien2013b():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/obrien2013_coccosphere_size.csv")
    #read obrien size data
    #apply synonyms
    d = rename_synonyms(d)
    d['sd'] = (d['max']-d['min'])/4
    d = d[['species', 'sd', 'mean']]
    d['reference'] = "obrien2013b"
    d['method'] = 'literature coccosphere'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/obrien2013b.csv", index=False)

def pre_devries2024():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/devries2024_volumes.csv")
    d = d[['species', 'mean', 'std']]
    d['sd'] = d['std']
    d['reference'] = 'devries2024'
    d['method'] = 'literature morphometrics'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/devries2024.csv", index=False)

def pre_sheward2024():

    #read sheward size data
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/sheward2024_volumes.csv")

    d = d[d['Spheres Y/N?']!="Flattened"]
    d['Species'] = d['Species'].str.strip()
    d['Genus'] = d['Genus'].str.strip()
    d['species'] = d['Genus'] + " " + d['Species'] 

    d = d[['species', 'Estimated cell volume', 'PIC pg C']]

    d = rename_synonyms(d, remove_duplicate=False)

    d = d.groupby(by="species").agg(["mean", "std"]).reset_index()
    d['sd'] = np.round(d['Estimated cell volume']['std'], 1)
    d['mean'] = np.round(d['Estimated cell volume']['mean'], 1)
    d = d[['species', 'sd', 'mean']]

    #rename because column names are tuples for some reason??
    d.columns = ['species', 'sd', 'mean']

    d['reference'] = 'sheward2024'
    d['method'] = "AMT morphometrics"
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/sheward2024.csv", index=False)


def pre_young2024():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/young2024.csv")
    d['mean'] = d['avg(volume)']
    d['sd'] = d['stddev(volume)']
    d['reference'] = 'young2024'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/young2024.csv", index=False)


def pre_gafar2019():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/gafar2019.csv")
    d['mean'] = d['cell volume']
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'gafar2019'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/gafar2019.csv", index=False)

def pre_fiorini2011():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/fiorini2011.csv")
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'fiorini2011'
    d['method'] = 'coulter counter'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/fiorini2011.csv", index=False)

def pre_gerecht2018():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/gerecht2018.csv")
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'gerecht2018'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/gerecht2018.csv", index=False)

def pre_mullin1966():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/mullin1966.csv")
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'mullin1966'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/mullin1966.csv", index=False)

def pre_oviedo2014():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/oviedo2014.csv")
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'oviedo2014'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/oviedo2014.csv", index=False)

def pre_supraha2015():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/supraha2015.csv")
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'supraha2015'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/supraha2015.csv", index=False)

def pre_verity1992():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sizes/verity1992.csv")
    d = rename_synonyms(d)
    d = resample_size_d(d)
    d['reference'] = 'verity1992'
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/verity1992.csv", index=False)


pre_villiot2021a()
pre_villiot2021b()
pre_obrien2013a()
# pre_obrien2013b() #now excluded from dataset
pre_devries2024()
pre_sheward2024()
pre_young2024()
pre_gafar2019()
pre_fiorini2011()
pre_gerecht2018()
pre_mullin1966()
pre_oviedo2014()
pre_supraha2015()
pre_verity1992()

'''
TO DO (with resampling):

pre_fiorini2011()
pre_gerecht2018()
pre_mullin1966()
pre_oviedo2014()
pre_supraha2015()
pre_verity1992()

'''
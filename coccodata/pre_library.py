import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 

def rename_synonyms(d, index='species', remove_duplicate=True):
    d['species'] = d['species'].str.strip()

    with open('/home/phyto/CoccoData/classification/synonyms.yml', 'r') as f:
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

def pre_villiot2021a():    
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/viliot2021_cell_diameters.csv")
    d = rename_synonyms(d, ['species', 'strain'])

    #resample ehux
    RCC1731 = np.random.normal(d[d["strain"]=="RCC1731"]['mean'], d[d["strain"]=="RCC1731"]['std'], 10000)
    RCC1128 = np.random.normal(d[d["strain"]=="RCC1228"]['mean'], d[d["strain"]=="RCC1228"]['std'], 10000)
    
    ehux = {'species': ["Emiliania huxleyi"],
                'mean': [np.mean([RCC1128, RCC1731])],
                'std': [np.std([RCC1128, RCC1731])],
                'strain': ['NA']}
    ehux = pd.DataFrame(ehux)
    d = d[d["species"]!="Emiliania huxleyi"]
    
    #resample leptoporus
    RCC1130 = np.random.normal(d[d["strain"]=="RCC1130"]['mean'], d[d["strain"]=="RCC1130"]['std'], 10000)
    RCC1135 = np.random.normal(d[d["strain"]=="RCC1135"]['mean'], d[d["strain"]=="RCC1135"]['std'], 10000)
    
    leptoporus = {'species': ["Coccolithus leptoporus"],
                'mean': [np.mean([RCC1130, RCC1135])],
                'sd': [np.std([RCC1130, RCC1135])],
                'strain':['NA']}
    leptoporus = pd.DataFrame(leptoporus)
    d = d[d["species"]!="Calcidiscus leptoporus"]
    d = pd.concat([d, ehux, leptoporus])
    d['mean'] = np.round(d['mean'])
    d['sd'] = np.round(d['sd'])
    d = d[['species', 'mean', 'sd']]
    d['reference'] = "villiot2021a"
    d['method'] = 'light microscopy'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/viliot2021a.csv", index=False)

def pre_villiot2021b():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/villiot2021_literature_morphometric_volume.csv")
    d = rename_synonyms(d)
    d = d.rename(columns={'cell volume': "mean"})
    d['mean'] = np.round(d['mean'])
    d['reference'] = "villiot2021b"
    d['method'] = 'literature morphometrics'
    d['sd'] = None
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/viliot2021b.csv", index=False)

def pre_obrien2013a():
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/obrien_cell_diameters.csv")
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
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/obrien2013_coccosphere_size.csv")
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
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/devries2024_volumes.csv")
    d = d[['species', 'mean', 'std']]
    d['sd'] = d['std']
    d['reference'] = 'devries2024'
    d['method'] = 'literature morphometrics'
    d = d[['species', 'mean', 'sd', 'method', 'reference']] 
    d.to_csv("/home/phyto/CoccoData/data/sizes/devries2024.csv", index=False)

def pre_sheward2024():

    #read sheward size data
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/sheward2024_volumes.csv")

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



pre_villiot2021a()
pre_villiot2021b()
pre_obrien2013a()
pre_obrien2013b()
pre_devries2024()
pre_sheward2024()



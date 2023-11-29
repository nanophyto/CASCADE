import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 
import glob, os

Study = namedtuple("Study", ['id', 'mean', 'sd', 'method', 'species', 'genera', 'family', 'phase', 'alternate_phase'])

def def_grouping():

    with open('/home/phyto/CoccoData/data/classification/phases.yml', 'r') as f:
        phases = load(f, Loader=Loader)

    with open('/home/phyto/CoccoData/data/classification/family.yml', 'r') as f:
        families = load(f, Loader=Loader)

    d = pd.DataFrame.from_dict(phases, orient='index')
    d = d.rename_axis("species").reset_index()
    d['genera'] = d['species'].str.split(" ").str[0]

    inverse = {}
    for k,v in families.items():
        for x in v:
            inverse.setdefault(x, []).append(k)

    df = pd.DataFrame.from_dict(inverse, orient='index')
    df = df.rename_axis("genera").reset_index()
    df = df.rename(columns={0: "family"})
    d = pd.merge(d, df, on='genera', how="outer")
    d = d.where(pd.notnull(d), None)

    library = list(d.itertuples(name='species', index=False))

    return(d)

groups = def_grouping()

species = groups.species
species = {x for x in species if x is not None}
species_list = list(species)
species_list.sort()

def import_data(path):

    all_files = glob.glob(os.path.join(path, "*.csv"))

    d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
    d = d.fillna(0)

    return(d)

d = import_data("/home/phyto/CoccoData/data/sizes/")

list_of_studies = d['reference'].unique()

def size_and_method(d, study, species):  

    d = d[(d['species'] == species) & (d['reference']==study)]

    keys_to_extract = ['mean', 'sd', 'method']
    extracted_data = {}

    for key in keys_to_extract:
        try:
            extracted_data[key] = d.get(key, None).item()
        except:
            extracted_data[key] = None
    return extracted_data['mean'], extracted_data['sd'], extracted_data['method']

def classification(groups, species):
    
    groups = groups[groups['species'] == species]
    keys_to_extract = ['genera', 'family', 'phase', 'alternate_phase']
    extracted_data = {}

    for key in keys_to_extract:
        try:
            extracted_data[key] = groups.get(key, None).item()
        except:
            extracted_data[key] = None

    return extracted_data['genera'], extracted_data['family'], extracted_data['phase'], extracted_data['alternate_phase']

def fill_namedtuple(groups, species_list):

    studies = []
    for id in list_of_studies:
        for species in species_list: 
            genera, family, phase, alternate_phase = classification(groups, species)
            studies.append(Study(id, *size_and_method(d, id, species), species, genera, family, phase, alternate_phase))

    return(studies)

library = fill_namedtuple(groups, species_list)

def export_yml(library, path):

    spp_list = []
    #sort species list alphabetically:


    for i in range(len(species_list)):

        name = species_list[i]
        species_library =  [t for t in library  if t.species == name]

        sizes = {study.id:study._asdict() for study in species_library }
        for id in list_of_studies:
            del sizes[id]['id']
            del sizes[id]['species']
            del sizes[id]['genera']
            del sizes[id]['family']
            del sizes[id]['phase']
            del sizes[id]['alternate_phase']
        
        #d = asdict(dc, dict_factory=lambda x: {k: v for (k, v) in x if v is not None})


        species =  {name: {
                'genera': species_library[0].genera,
                'family': species_library[0].family,
               'phase': species_library[0].phase,
                'alternate_phase': species_library[0].alternate_phase,
                'size' : sizes
            }}
        spp_list.append(species)

    with open(path, 'w') as outfile:
        yaml.dump(spp_list, outfile, default_flow_style=False)

    print("exported yml to: " + str(path))

export_yml(library, '/home/phyto/CoccoData/library.yml')



#quick and dirty export (might be easier to read?)
# def export_yml(library, path):

#     spp_list = []

#     for i in range(len(library)):
#         species = library[i]._asdict()
#         spp_list.append(species)

#     with open(path, 'w') as outfile:
#         yaml.dump(spp_list, outfile, default_flow_style=False)

#     print("exported yml to: " + str(path))


# export_yml(library, '/home/phyto/CoccoData/library.yml')



def find_undefined_spp(library):
    for i in range(len(library)):
        ntpl = library[i]

        if (
            (ntpl.size.obrien2013a.mean ==None) and 
            (ntpl.size.obrien2013b.mean ==None) and 
            (ntpl.size.villiot2021a.mean==None) and 
            (ntpl.size.villiot2021b.mean==None) and 
            (ntpl.size.devries2024.mean==None) and 
            (ntpl.size.sheward2024.mean==None) ):
                print(str(ntpl.species))


find_undefined_spp(library)

print("fin")
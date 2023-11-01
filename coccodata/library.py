import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 

Study = namedtuple("Study", ['id', 'mean', 'sd', 'method', 'species'])



def def_grouping():

    with open('/home/phyto/CoccoData/classification/phases.yml', 'r') as f:
        phases = load(f, Loader=Loader)

    with open('/home/phyto/CoccoData/classification/family.yml', 'r') as f:
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

species_list = groups.species.tolist()



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

def import_villiot2021a():    
    d = pd.read_csv("/home/phyto/CoccoData/sizes/viliot2021_cell_diameters.csv")
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
                'std': [np.std([RCC1130, RCC1135])],
                'strain':['NA']}
    leptoporus = pd.DataFrame(leptoporus)

    d = d[d["species"]!="Calcidiscus leptoporus"]

    d = pd.concat([d, ehux, leptoporus])
    d['mean'] = np.round(d['mean'])
    d['std'] = np.round(d['std'])
    villiot2021 = d[['species', 'mean', 'std']]

    villiot2021['reference'] = "villiot2021a"
    villiot2021['method'] = 'light microscopy'

    return(villiot2021)

villiot2021a = import_villiot2021a()




def import_villiot2021b():
    d = pd.read_csv("/home/phyto/CoccoData/sizes/villiot2021_literature_morphometric_volume.csv")
    d = rename_synonyms(d)
    d = d.rename(columns={'cell volume': "mean"})
    d['mean'] = np.round(d['mean'])

    d['reference'] = "villiot2021b"
    d['method'] = 'literature morphometrics'

    return(d)

villiot2021b = import_villiot2021b()



def import_obrien2013a():
    d = pd.read_csv("/home/phyto/CoccoData/sizes/obrien_cell_diameters.csv")
    
    #read obrien size data
    #apply synonyms
    d = rename_synonyms(d)

    d = d.rename(columns={'diameter': "mean"})
    d['mean'] = (1/6)*math.pi*(d['mean']**3)
    d['mean'] = np.round(d['mean'])
    d['reference'] = 'obrien2013a'
    d['std'] = None
    d['method'] = 'light microscopy'
    return(d)

obrien2013a = import_obrien2013a()



def import_obrien2013b():
    d = pd.read_csv("/home/phyto/CoccoData/sizes/obrien2013_coccosphere_size.csv")
    #read obrien size data
    #apply synonyms
    d = rename_synonyms(d)

    d['std'] = (d['max']-d['min'])/4
    d = d[['species', 'std', 'mean']]
    d['reference'] = "obrien2013b"
    d['method'] = 'literature coccosphere'
    return(d)

obrien2013b = import_obrien2013b()



def import_devries2024():
    d = pd.read_csv("/home/phyto/CoccoData/sizes/devries2024_volumes.csv")
    d = d[['species', 'mean', 'std']]
    d['reference'] = 'devries2024'
    d['method'] = 'literature morphometrics'

    return(d)

devries2024 = import_devries2024()

def import_sheward2024():

    #read sheward size data
    d = pd.read_csv("/home/phyto/CoccoData/sizes/sheward2024_volumes.csv")

    d = d[d['Spheres Y/N?']!="Flattened"]
    d['Species'] = d['Species'].str.strip()
    d['Genus'] = d['Genus'].str.strip()
    d['species'] = d['Genus'] + " " + d['Species'] 

    d = d[['species', 'Estimated cell volume', 'PIC pg C']]

    d = rename_synonyms(d, remove_duplicate=False)


    d = d.groupby(by="species").agg(["mean", "std"]).reset_index()
    d['std'] = np.round(d['Estimated cell volume']['std'], 1)
    d['mean'] = np.round(d['Estimated cell volume']['mean'], 1)
    sheward2024 = d[['species', 'std', 'mean']]

    #rename because column names are tuples for some reason??
    sheward2024.columns = ['species', 'std', 'mean']

    sheward2024['reference'] = 'sheward2024'
    sheward2024['method'] = "AMT morphometrics"

    return(sheward2024)

sheward2024 = import_sheward2024()


d = pd.concat([obrien2013a, obrien2013b, sheward2024,
                villiot2021a, villiot2021b, devries2024])


def def_sizes(d, study, species):
    d = d[(d['species'] == species) & (d['reference']==study)]
    try:
        method = d.get('method', None).item()
    except:
        method = None
    try:
        d_mean = d.get('mean', None).item()
    except:
        d_mean = None
    try:
        d_std = d.get('std', None).item()
    except:
        d_std = None
    return d_mean, d_std, method


def fill_namedtuple(groups, species_list):

    studies = []
    for id in ['obrien2013a', 'obrien2013b', 'sheward2024',
                'villiot2021a', 'villiot2021b', 'devries2024']:
        for species in species_list: 
            studies.append(Study(id, *def_sizes(d, id, species), species))

    return(studies)

library = fill_namedtuple(groups, species_list)


# def export_yml(library, path):

#     spp_list = []

#     for i in range(len(species_list)):
#         ntpl = library[i]
#         name = ntpl.species
#         species =  {name: {
# #                'genera': ntpl.genera,
# #                'family': ntpl.family,
# #                'phase': ntpl.phase,
# #                'alternate_phase': ntpl.alternate_phase,
#                 'size' : {study.id:study._asdict() for study in species_list}
#             }}
#         spp_list.append(species)

#     with open(path, 'w') as outfile:
#         yaml.dump(spp_list, outfile, default_flow_style=False)

#     print("exported yml to: " + str(path))

# export_yml(library, '/home/phyto/CoccoData/library.yml')




def export_yml(library, path):

    spp_list = []

    for i in range(len(library)):
        species = library[i]._asdict()
        spp_list.append(species)

    with open(path, 'w') as outfile:
        yaml.dump(spp_list, outfile, default_flow_style=False)

    print("exported yml to: " + str(path))


export_yml(library, '/home/phyto/CoccoData/library.yml')



def find_undefined_spp(library):
    for i in range(len(library)):
        ntpl = library[i]

        if (
            (convert_float(ntpl.size.obrien2013a.mean) ==None) and 
            (convert_float(ntpl.size.obrien2013b.mean) ==None) and 
            (convert_float(ntpl.size.villiot2021a.mean)==None) and 
            (convert_float(ntpl.size.villiot2021b.mean)==None) and 
            (convert_float(ntpl.size.devries2024.mean)==None) and 
            (convert_float(ntpl.size.sheward2024.mean)==None) ):
                print(str(ntpl.species))


find_undefined_spp(library)

print("fin")
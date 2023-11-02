import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 

Study = namedtuple("Study", ['id', 'mean', 'sd', 'method', 'species', 'genera', 'family', 'phase', 'alternate_phase'])



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


def size_and_method(d, study, species):  

    d = d[(d['species'] == species) & (d['reference']==study)]

    keys_to_extract = ['mean', 'std', 'method']
    extracted_data = {}

    for key in keys_to_extract:
        try:
            extracted_data[key] = d.get(key, None).item()
        except:
            extracted_data[key] = None

    return extracted_data['mean'], extracted_data['std'], extracted_data['method']

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
    for id in ['obrien2013a', 'obrien2013b', 'sheward2024',
                'villiot2021a', 'villiot2021b', 'devries2024']:
        for species in species_list: 
            genera, family, phase, alternate_phase = classification(groups, species)

            studies.append(Study(id, *size_and_method(d, id, species), species, genera, family, phase, alternate_phase))

    return(studies)

library = fill_namedtuple(groups, species_list)

def export_yml(library, path):

    spp_list = []

    for i in range(len(species_list)):

        name = species_list[i]
        species_library =  [t for t in library  if t.species == name]

        sizes = {study.id:study._asdict() for study in species_library }

        for id in ['obrien2013a', 'obrien2013b', 'sheward2024',
                        'villiot2021a', 'villiot2021b', 'devries2024']:
            del sizes[id]['id']
            del sizes[id]['species']
            del sizes[id]['genera']
            del sizes[id]['family']
            del sizes[id]['phase']
            del sizes[id]['alternate_phase']
            
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
            (convert_float(ntpl.size.obrien2013a.mean) ==None) and 
            (convert_float(ntpl.size.obrien2013b.mean) ==None) and 
            (convert_float(ntpl.size.villiot2021a.mean)==None) and 
            (convert_float(ntpl.size.villiot2021b.mean)==None) and 
            (convert_float(ntpl.size.devries2024.mean)==None) and 
            (convert_float(ntpl.size.sheward2024.mean)==None) ):
                print(str(ntpl.species))


find_undefined_spp(library)

print("fin")
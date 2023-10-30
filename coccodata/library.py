import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math 

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
    print(df)
    d = pd.merge(d, df, on='genera', how="outer")
    d = d.where(pd.notnull(d), None)

    library = list(d.itertuples(name='species', index=False))

    return(d)

groups = def_grouping()

species_list = groups.species



def rename_synonyms(d, index='species'):
    d['species'] = d['species'].str.strip()

    with open('/home/phyto/CoccoData/classification/synonyms.yml', 'r') as f:
        groupings = load(f, Loader=Loader)

    species = d['species']
    synonym_dict = {species:k
        for k, v in groupings.items()
        for species in v['alt']}

    d = d.replace(synonym_dict)

    #take mean for duplicate entries
    d = d.groupby(index).mean().reset_index()    

    return(d)

def import_villiot2021():    
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

    return(villiot2021)

villiot2021 = import_villiot2021()




def import_villiot2021b():
    d = pd.read_csv("/home/phyto/CoccoData/sizes/villiot2021_literature_morphometric_volume.csv")

    #read obrien size data
    d = rename_synonyms(d)
    d = d.rename(columns={'cell volume': "mean"})
    d['mean'] = np.round(d['mean'])
    return(d)

villiot2021b = import_villiot2021b()



def import_obrien2013():
    d = pd.read_csv("/home/phyto/CoccoData/sizes/obrien_cell_diameters.csv")
    
    #read obrien size data
    #apply synonyms
    d = rename_synonyms(d)

    d = d.rename(columns={'diameter': "mean"})
    d['mean'] = 1/6*math.pi*(d['mean']**3)
    d['mean'] = np.round(d['mean'])
    return(d)

obrien2013 = import_obrien2013()


def import_sheward2024():

    #read sheward size data
    d = pd.read_csv("/home/phyto/CoccoData/sizes/sheward2024_volumes.csv")
    d['Species'] = d['Species'].str.strip()
    d['Genus'] = d['Genus'].str.strip()
    d['species'] = d['Genus'] + " " + d['Species'] 

    d = d[['species', 'Estimated cell volume', 'PIC pg C']]

    d = rename_synonyms(d)


    d_mean = d.groupby(by="species").mean()

    d_std = d.groupby(by="species").std()

    d_mean = d_mean.rename(columns={'PIC pg C': "PIC_mean",
                                            'Estimated cell volume': "mean"}).reset_index()

    d_std = d_std.rename(columns={'PIC pg C': "PIC_std",
                                            'Estimated cell volume': "std"})

    d_std['std'] = np.round(d_std['std'])
    d_mean['mean'] = np.round(d_mean['mean'])
    d = pd.concat([d_mean[['species', 'mean']], d_std['std'] ], axis=1)

    return(d)
sheward2024 = import_sheward2024()


                                            
def def_sizes(d, species):

    d = d[d['species']==species]

    print(d)

    try:
        d_mean = d['mean']
    except:
        d_mean = None
    try:
        d_std = d['std']
    except:
        d_std = None

    return(d_mean, d_std)


print(groups['species'])


def create_namedtuple(groups, species):

    obrien2013_values = namedtuple('obrien2013', 'mean std method')
    sheward2024_values = namedtuple('sheward2024', 'mean std method')
    villiot2021a_values = namedtuple('villiot2021a', 'mean std method')
    villiot2021b_values = namedtuple('villiot2021b', 'mean std method')

    sizes = namedtuple('size', ['obrien2013', 'sheward2024',
                                'villiot2021a', 'villiot2021b'])

    library = namedtuple('library', ['species', 'genera', 
            'family', 'phase', 'alternate_phase', 'size'])

    groups = groups[groups['species']==species]

    obrien2013_mean, obrien2013_std  = def_sizes(obrien2013, species)
    sheward2024_mean, sheward2024_std  = def_sizes(sheward2024, species)
    villiot2021a_mean, villiot2021a_std  = def_sizes(villiot2021, species)
    villiot2021b_mean, villiot2021b_std = def_sizes(villiot2021b, species)

    ntpl = library(
        groups['species'], 
        groups['genera'], 
        groups['family'],  
        groups['phase'], 
        groups['alternate_phase'],
        sizes(
            obrien2013_values(
                obrien2013_mean, 
                obrien2013_std,
                "light microscopy"
            ),
            sheward2024_values(
                sheward2024_mean, 
                sheward2024_std,
                "AMT morphometrics"
            ),
            villiot2021a_values(
                villiot2021a_mean, 
                villiot2021a_std,
                "light microscopy"
            ),
            villiot2021b_values(
                villiot2021b_mean, 
                villiot2021b_std,
                "literature morphometrics"
            ),

        )
    )

    return(ntpl)

library = []

for i in range(len(species_list)):
    library.append(create_namedtuple(groups, species_list[i]))

def convert_float(a):
    try:
        b = float(a)
    except:
        b = None
    return(b)

def convert_str(a):
    try:
        b = str(a.iloc[0])
    except:
        b = None
    return(b)

def export_yml(library, path):

    spp_list = []

    for i in range(len(library)):
        ntpl = library[i]
        print(ntpl.species)
        name = convert_str(ntpl.species)
        print(name)
        species =  {name: {
                'genera': convert_str(ntpl.genera),
                'family': convert_str(ntpl.family),
                'phase': convert_str(ntpl.phase),
                'alternate_phase': convert_str(ntpl.alternate_phase),
                'size':{
                    'obrien2013':{
                        'mean': convert_float(ntpl.size.obrien2013.mean),
                        'std': convert_float(ntpl.size.obrien2013.std),
                        'method': ntpl.size.obrien2013.method
                        },
                    'villiot2021a':{
                        'mean': convert_float(ntpl.size.villiot2021a.mean),
                        'std': convert_float(ntpl.size.villiot2021a.std),
                        'method': ntpl.size.villiot2021a.method
                        },
                    'villiot2021b':{
                        'mean': convert_float(ntpl.size.villiot2021b.mean),
                        'std': convert_float(ntpl.size.villiot2021b.std),
                        'method': ntpl.size.villiot2021b.method
                        },                
                    'sheward2024':{
                        'mean': convert_float(ntpl.size.sheward2024.mean),
                        'std': convert_float(ntpl.size.sheward2024.std),
                        'method': ntpl.size.sheward2024.method
                        }
                    }
                }
            }
        spp_list.append(species)

    with open(path, 'w') as outfile:
        yaml.dump(spp_list, outfile, default_flow_style=False)

    print("exported yml to: " + str(path))


export_yml(library, '/home/phyto/CoccoData/library.yml')



print("fin")



# #check which HET species have undefined volumes:
# volume_check = d[d['phase']=="HET"]
# volume_check['sum'] = volume_check[['size_mean_sheward2024', 
#                                     'size_mean_obrien2013']].mean(axis=1)
# test = volume_check[['sum', 'species']] 
# print(test[test.isna().any(axis=1)]['species'])


# #check which HOL species have undefined volumes:
# volume_check = d[d['phase']=="HOL"]
# volume_check = volume_check[volume_check['2N']=="unknown"]

# volume_check['sum'] = volume_check[['size_mean_sheward2024', 
#                                     'size_mean_obrien2013']].mean(axis=1)
# test = volume_check[['sum', 'species']] 
# print(test[test.isna().any(axis=1)]['species'])


# d.to_csv("/home/phyto/CoccoData/library.csv", index=False)

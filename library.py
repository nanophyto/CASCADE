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


def import_villiot2021():    
    villiot_size = pd.read_csv("/home/phyto/CoccoData/sizes/viliot2021_cell_diameters.csv")

    #resample ehux
    RCC1731 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1731"]['mean'], villiot_size[villiot_size["strain"]=="RCC1731"]['std'], 10000)
    RCC1128 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1228"]['mean'], villiot_size[villiot_size["strain"]=="RCC1228"]['std'], 10000)
    
    ehux = {'species': ["Emiliania huxleyi"],
                'mean': [np.mean([RCC1128, RCC1731])],
                'std': [np.std([RCC1128, RCC1731])],
                'strain': ['NA']}
    ehux = pd.DataFrame(ehux)
    villiot_size = villiot_size[villiot_size["species"]!="Emiliania huxleyi"]
    
    #resample leptoporus
    RCC1130 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1130"]['mean'], villiot_size[villiot_size["strain"]=="RCC1130"]['std'], 10000)
    RCC1135 = np.random.normal(villiot_size[villiot_size["strain"]=="RCC1135"]['mean'], villiot_size[villiot_size["strain"]=="RCC1135"]['std'], 10000)
    
    leptoporus = {'species': ["Coccolithus leptoporus"],
                'mean': [np.mean([RCC1130, RCC1135])],
                'std': [np.std([RCC1130, RCC1135])],
                'strain':['NA']}
    leptoporus = pd.DataFrame(leptoporus)

    villiot_size = villiot_size[villiot_size["species"]!="Calcidiscus leptoporus"]

    villiot_size = pd.concat([villiot_size, ehux, leptoporus])

    villiot2021 = villiot_size[['species', 'mean', 'std']]

    return(villiot2021)

villiot2021 = import_villiot2021()




def import_obrien2013():
    
    #read obrien size data
    obrien_size = pd.read_csv("/home/phyto/CoccoData/sizes/obrien_cell_diameters.csv")
    #apply synonyms
    obrien_size['species'] = obrien_size['species'].str.strip()

    with open('/home/phyto/CoccoData/classification/synonyms.yml', 'r') as f:
        groupings = load(f, Loader=Loader)

    species = obrien_size['species']
    dict = {species:k
        for k, v in groupings.items()
        for species in v['alt']}

    obrien_size = obrien_size.replace(dict)

    df = obrien_size.rename(columns={'diameter': "mean"})

    #take mean for duplicate entries
    df = df.groupby('species').mean().reset_index()    

    ##########
    # STILL NEED TO ESTIMATE VOLUME FROM DIAMETER
    ##########

    #merge size data
    return(df)

obrien2013 = import_obrien2013()


def import_sheward2024():

    #read sheward size data
    sheward2024 = pd.read_csv("/home/phyto/CoccoData/sheward2024.csv")
    sheward2024['Species'] = sheward2024['Species'].str.strip()
    sheward2024['Genus'] = sheward2024['Genus'].str.strip()
    sheward2024['species'] = sheward2024['Genus'] + " " + sheward2024['Species'] 

    species = sheward2024['species']
    dict = {species:k
        for k, v in groupings.items()
        for species in v['alt']}

    sheward2024 = sheward2024.replace(dict)

    sheward2024 = sheward2024[['species', 'Estimated cell volume', 'PIC pg C']]

    sheward2024_mean = sheward2024.groupby(by="species").mean()

    sheward2024_std = sheward2024.groupby(by="species").std()

    sheward2024_mean = sheward2024_mean.rename(columns={'PIC pg C': "PIC_mean",
                                            'Estimated cell volume': "size_mean_sheward2024"})

    sheward2024_std = sheward2024_std.rename(columns={'PIC pg C': "PIC_std",
                                            'Estimated cell volume': "size_std_sheward2024"})



                                            
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


def create_namedtuple(groups, villiot2021, obrien2013, species):

    mean_values = namedtuple('mean', 'villiot2021 obrien2013 sheward2024')
    std_values = namedtuple('std', 'villiot2021 obrien2013 sheward2024')

    sizes = namedtuple('size', ['mean', 'std'])

    library = namedtuple('library', ['species', 'genera', 
            'family', 'phase', 'alternate_phase', 'size'])

    groups = groups[groups['species']==species]

    villiot_mean, villiot_std  = def_sizes(villiot2021, species)
    obrien_mean, obrien_std  = def_sizes(obrien2013, species)
    sheward_mean, sheward_std  = def_sizes(obrien2013, species)

    ntpl = library(
        groups['species'], 
        groups['genera'], 
        groups['family'],  
        groups['phase'], 
        groups['alternate_phase'],
        sizes(
            mean_values(
                villiot_mean, 
                obrien_mean, 
                sheward_mean), 
            std_values(
                villiot_std, 
                obrien_std, 
                sheward_std)
            )
        )

    return(ntpl)

library = []

for i in range(len(species_list)):
    library.append(create_namedtuple(groups, villiot2021, obrien2013, species_list[i]))

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
        name = convert_str(ntpl.species)
        print(name)
        species =  {name: {
                'genera': convert_str(ntpl.genera),
                'family': convert_str(ntpl.family),
                'phase': convert_str(ntpl.phase),
                'alternate_phase': convert_str(ntpl.alternate_phase),
                'size':{
                    'mean':{
                        'obrien': convert_float(ntpl.size.mean.obrien2013),
                        'villiot':convert_float(ntpl.size.mean.villiot2021),
                        'sheward':convert_float(ntpl.size.mean.sheward2024)
                    },
                    'std':{
                        'obrien':convert_float(ntpl.size.std.obrien2013),
                        'villiot':convert_float(ntpl.size.std.villiot2021),
                        'sheward':convert_float(ntpl.size.std.sheward2024)
                    }
                }
            }
        }
        spp_list.append(species)

    with open(path, 'w') as outfile:
        yaml.dump(spp_list, outfile)

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

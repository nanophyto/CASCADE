import numpy as np 
import pandas as pd
import math

def volume_cone(d, h):
    volume = math.pi*((d/2)**2)*(h/3)
    return(volume)

def volume_sphere(d):
    volume = (1/6)*math.pi*(d**3)
    return(volume)


sample_n = 10000

def aurisinae():
    coccosphere_length_min = 23
    coccosphere_length_max = 26

    coccosphere_width_min = 13
    coccosphere_width_max = 17

    coccolith_width_min = 1.1
    coccolith_width_max = 1.2


    coccosphere_length_sd = (coccosphere_length_max-coccosphere_length_min)/4
    coccosphere_length_mean = (coccosphere_length_max+coccosphere_length_min)/2

    coccosphere_width_sd = (coccosphere_width_max-coccosphere_width_min)/4
    coccosphere_width_mean = (coccosphere_width_max+coccosphere_width_min)/2

    coccolith_thickness_sd = (coccolith_width_max-coccolith_width_min)/4
    coccolith_thickness_mean = (coccolith_width_max+coccolith_width_min)/2


    cell_length = np.clip(np.random.normal(coccosphere_length_mean, coccosphere_length_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)

    cell_width = np.clip(np.random.normal(coccosphere_width_mean, coccosphere_width_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)


    cell_volume = volume_cone(cell_width, cell_length)

    cell_mean = np.mean(cell_volume)
    cell_std = np.std(cell_volume)


    d = {'species':['Syracosphaera aurisinae'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['cone'],
        'ref': ["lecal-schlauder1951"]}

    d = pd.DataFrame(d)
    return(d)


def caudatus():
    coccosphere_length_min = 26
    coccosphere_length_max = 36

    coccosphere_width_min = 3.5
    coccosphere_width_max = 4

    coccolith_thickness = 1


    coccosphere_length_sd = (coccosphere_length_max-coccosphere_length_min)/4
    coccosphere_length_mean = (coccosphere_length_max+coccosphere_length_min)/2

    coccosphere_width_sd = (coccosphere_width_max-coccosphere_width_min)/4
    coccosphere_width_mean = (coccosphere_width_max+coccosphere_width_min)/2


    cell_length = np.clip(np.random.normal(coccosphere_length_mean, coccosphere_length_sd, sample_n), 0, None) - 2*coccolith_thickness

    cell_width = np.clip(np.random.normal(coccosphere_width_mean, coccosphere_width_sd, sample_n), 0, None) - 2*coccolith_thickness


    cell_volume = volume_cone(cell_width, cell_length)

    cell_mean = np.mean(cell_volume)
    cell_std = np.std(cell_volume)


    d = {'species':['Calciopappus caudatus'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['cone'],
        'ref': ["gaarder1971"]}

    d = pd.DataFrame(d)
    return(d)




def borealis():
    coccosphere_min = 6.5 
    coccosphere_max = 8.2
    coccolith_thickness = 1

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Syracosphaera borealis'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['okada1977']}

    d = pd.DataFrame(d)

    return(d)



def adenensis():
    coccosphere_min = 5.5 
    coccosphere_max = 8.5
    coccolith_thickness = 1

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Sphaerocalyptra adenensis'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)



def blokii():
    coccosphere = 6
    coccolith_thickness = 1

    diameters = coccosphere - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)

    d = {'species':['Calicasphaera blokii'], 
        'mean': [cell_mean],
        'std': 0,
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)



def concava():
    coccosphere = 6
    coccolith_thickness = 1.3

    diameters = coccosphere - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)

    d = {'species':['Calicasphaera concava'], 
        'mean': [cell_mean],
        'std': 0,
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)



def gaudii_POL():
    coccosphere_min = 5.6 
    coccosphere_max = 10.6
    coccolith_thickness = 1

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Alisphaera gaudii POL'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)



def pienaarii():
    coccosphere_min = 6 
    coccosphere_max = 11
    coccolith_thickness = 0.9

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Helladosphaera pienaarii'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['norris1985, kleijne1991']}

    d = pd.DataFrame(d)

    return(d)


def strigilis():
    coccosphere_min = 6 
    coccosphere_max = 9
    coccolith_thickness = 1

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*coccolith_thickness

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Syracosphaera strigilis'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)



def formosus():
    coccosphere_min = 4.5 
    coccosphere_max = 7.5
    coccolith_thickness_min = 0.7
    coccolith_thickness_max = 0.8

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2
    coccolith_thickness_sd = (coccolith_thickness_max-coccolith_thickness_min)/4
    coccolith_thickness_mean = (coccolith_thickness_max+coccolith_thickness_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)

    volumes = volume_sphere(diameters)

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Ophiaster formosus'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)

def type_5():

    coccosphere_min = 6
    coccosphere_max = 7

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters =  np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None)

    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Pappomonas sp. type 5'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)
    return(d)



def squamosa():
    coccosphere = 6 #um
    coccolith_thickness_min = 0.9
    coccolith_thickness_max = 1.3

    coccolith_thickness_sd = (coccolith_thickness_max-coccolith_thickness_min)/4
    coccolith_thickness_mean = (coccolith_thickness_max+coccolith_thickness_min)/2

    diameters = 5-(2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None))

    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Syracosphaera squamosa'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['kleijne2009']}

    d = pd.DataFrame(d)

    return(d)

def reniformis():
    coccosphere_min = 6 
    coccosphere_max = 8
    coccolith_thickness_min = 1
    coccolith_thickness_max = 1.2

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2
    coccolith_thickness_sd = (coccolith_thickness_max-coccolith_thickness_min)/4
    coccolith_thickness_mean = (coccolith_thickness_max+coccolith_thickness_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)


    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Syracosphaera reniformis'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['klijne2009']}

    d = pd.DataFrame(d)

    return(d)

def sphaeroidea_hol():

    coccosphere_min = 5.5 
    coccosphere_max = 12

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None)

    volumes = volume_sphere(diameters)*0.9

    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Calyptrosphaera sphaeroidea HOL'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)


def arctica():


    cell_mean = volume_sphere(5)

    d = {'species':['Wigwamma antarctica'], 
        'mean': [cell_mean],
        'std': [None],
        'shape':['sphere'],
        'ref': ['thomsen1988']}

    d = pd.DataFrame(d)
    return(d)

def pienaarii():

    coccosphere_min = 6 
    coccosphere_max = 11
    coccolith_thickness_min = 0.7
    coccolith_thickness_max = 0.9

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2
    coccolith_thickness_sd = (coccolith_thickness_max-coccolith_thickness_min)/4
    coccolith_thickness_mean = (coccolith_thickness_max+coccolith_thickness_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)


    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Helladosphaera pienaarii'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['norris1985']}

    d = pd.DataFrame(d)

    return(d)


def cristatus():

    coccosphere_min =  7 #obrien2013 and refs within
    coccosphere_max = 18.9 #obrien2013 and refs within

    coccolith_thickness_min = 0.2 #young1998
    coccolith_thickness_max = 0.4 #young1998

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2
    coccolith_thickness_sd = (coccolith_thickness_max-coccolith_thickness_min)/4
    coccolith_thickness_mean = (coccolith_thickness_max+coccolith_thickness_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)


    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Ceratolithus cristatus'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['obrien2013 and young1998']}

    d = pd.DataFrame(d)

    return(d)




def hydroideus():

    coccolith_thickness_min = 0.7 
    coccolith_thickness_max = 0.9 

    coccosphere_sd = 0
    coccosphere_mean = 6
    coccolith_thickness_sd = (coccolith_thickness_max-coccolith_thickness_min)/4
    coccolith_thickness_mean = (coccolith_thickness_max+coccolith_thickness_min)/2

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)


    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Ophiaster hydroideus'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['cros2002']}

    d = pd.DataFrame(d)

    return(d)

def marsilii():
 
    coccosphere_min = 6.5  
    coccosphere_max = 8.5 

    coccosphere_sd = (coccosphere_max-coccosphere_min)/4
    coccosphere_mean = (coccosphere_max+coccosphere_min)/2
    coccolith_thickness_sd = 0
    coccolith_thickness_mean = 0.3

    diameters = np.clip(np.random.normal(coccosphere_mean, coccosphere_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)


    volumes = volume_sphere(diameters)
    cell_mean = np.round(np.mean(volumes), 1)
    cell_std = np.round(np.std(volumes), 1)

    d = {'species':['Zygosphaera marsilii'], 
        'mean': [cell_mean],
        'std': [cell_std],
        'shape':['sphere'],
        'ref': ['borsetti1976 and cros2002' ]}

    d = pd.DataFrame(d)

    return(d)   

d_reniformis = reniformis()
d_squamosa = squamosa()
d_formosus = formosus()
d_type5 = type_5()
d_aurisinae =  aurisinae()
d_arctica = arctica()
d_pienaarii = pienaarii()
d_sphaeroidea_hol = sphaeroidea_hol()
d_cristatus = cristatus()
d_hydroideus = hydroideus()
d_marsilii = marsilii()


d_caudatus = caudatus()
d_borealis = borealis()
d_adenensis = adenensis()
d_blokii = blokii()
d_concava = concava()
d_gaudii_POL = gaudii_POL()
d_pienaarii = pienaarii()
d_strigilis = strigilis()


d = pd.concat([d_type5, d_aurisinae, d_sphaeroidea_hol,
    d_reniformis, d_squamosa, d_arctica, d_pienaarii, d_cristatus, d_hydroideus,
    d_marsilii, d_caudatus, d_borealis, d_adenensis, d_blokii, d_concava, d_gaudii_POL,
    d_pienaarii, d_strigilis])

#read observation counts:

counts = pd.read_csv("/home/phyto/CoccoData/data/species_list_full.csv")
d.reset_index(inplace=True, drop=True)

d = pd.merge(d, counts, on="species")

d = d[['species','mean', 'std', 'shape', 'ref',  'count']]
d['std'] = d['std'].astype('float')

d.rename(columns={"species": "Species", "mean": "Mean ESD", "std": "SD ESD", 
                  "shape":"Cell shape", "ref": "Reference", "count":"Abundance obs."}, inplace=True)


print(d.to_latex(index=False,
                  formatters={"name": str.upper},
                  float_format="{:.1f}".format,
)) 


d['reference'] = 'devries2025'
d.to_csv("/home/phyto/CoccoData/data/sizes/devries2024.csv", index=False)

print("fin")



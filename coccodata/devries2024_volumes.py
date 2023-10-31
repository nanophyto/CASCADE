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

    d = {'species':['Pappomonas type 5'], 
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
        'ref': ['klijne2009']}

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



d_reniformis = reniformis()
d_squamosa = squamosa()
d_formosus = formosus()
d_type5 = type_5()
d_aurisinae =  aurisinae()
d_arctica = arctica()
d_pienaarii = pienaarii()
d_sphaeroidea_hol = sphaeroidea_hol()

d = pd.concat([d_formosus, d_type5, d_aurisinae, d_sphaeroidea_hol,
    d_reniformis, d_squamosa, d_arctica, d_pienaarii ])

d.to_csv("/home/phyto/CoccoData/sizes/devries2024_volumes.csv", index=False)

print("fin")



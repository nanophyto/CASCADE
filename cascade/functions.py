import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
from joblib import Parallel, delayed
np.random.seed(2)
import math

def rename_synonyms(d, classification_path='../data/classification/synonyms.yml',
                    verbose=0, index='species', 
                    remove_duplicate=True, take_sum = False, 
                    check_synonyms = True):
    d['species'] = d['species'].str.strip()

    with open(classification_path, 'r') as f:
        groupings = load(f, Loader=Loader)

    synonym_dict = {species:k
        for k, v in groupings.items()
        for species in v}

    d = d.replace(synonym_dict)

    #take mean for duplicate entries
    if remove_duplicate == True: 
        if take_sum == False:
            d = d.groupby(index).mean().reset_index()    
        else:
            d = d.groupby(index).sum().reset_index()    


    try:
        d = d[d['species']!='Thoracosphaera heimii']
    except:
        None

    try:
        d = d[d['species']!='Reticulofenestra sessilis']
    except:
        None

    #check if final species are in species library
    species_library = {v: k for k, v in synonym_dict.items()}

    species_observed = d['species']

    if check_synonyms == True:
        if (set(species_observed).issubset(species_library)):
            if verbose==1:
                print("all species are defined")
        else:
            raise ValueError("undefined species:" + str(set(species_observed).difference(species_library)))
    #        print("species missing: ")
    return(d)


def bayes_bootstrap(df, K=1000):

    def bayes_boot(df, estimator, seed=1):
        np.random.seed(seed)
        w = np.random.dirichlet(np.ones(len(df))*4, 1)[0]
        result = estimator(df, weights=w)
        return result

    r = Parallel(n_jobs=8)(delayed(bayes_boot)(df, np.average, seed=i) for i in range(K))

    return r


def ratio_bootstrap(x, y, K=100):

    def estimator(x, y):
        z = x/y
        return(z)

    def boot(x, y, seed=1):

        a = np.random.choice(x, 1)
        b = np.random.choice(y, 1)

        result = estimator(a, b)
        return result

    r = Parallel(n_jobs=8)(delayed(boot)(x, y, seed=i) for i in range(K))

    return np.concatenate(r).ravel()


def diameter_to_volume(d):
    v = (1/6) * np.pi * d**3
    return(v)



def abundance_refs_table(path, tex=False):
    d = pd.read_csv(path)

    refs = d['Reference'].unique()

    table = []

    for i in range (0, len(refs)):
        df = d[d['Reference']==refs[i]]
        t = pd.DataFrame({
            'Reference':[refs[i]],
            'Method':df['Method'].iloc[0],
            'Samples': len(df),
            'Survey Period': str(np.min(df['Year'])) + "-" + str(np.max(df['Year']))
        })
        table.append(t)
    
    table = pd.concat(table).sort_values(by=['Reference'])

    if tex:
        print(table.to_latex(index=False,
                        formatters={"name": str.upper},
                        float_format="{:.1f}".format,
        )) 


    return(table)


def add_months_since_solstice(d, month="Month"):
    d['months since winter solstice'] = d[month]
    d.loc[(d['Latitude'] <0) & (d['Month'] == 1),'months since winter solstice']= 7
    d.loc[(d['Latitude'] <0) & (d['Month'] == 2),'months since winter solstice']= 8
    d.loc[(d['Latitude'] <0) & (d['Month'] == 3),'months since winter solstice']= 9
    d.loc[(d['Latitude'] <0) & (d['Month'] == 4),'months since winter solstice']= 10
    d.loc[(d['Latitude'] <0) & (d['Month'] == 5),'months since winter solstice']= 11
    d.loc[(d['Latitude'] <0) & (d['Month'] == 6),'months since winter solstice']= 12
    d.loc[(d['Latitude'] <0) & (d['Month'] == 7),'months since winter solstice']= 1
    d.loc[(d['Latitude'] <0) & (d['Month'] == 8),'months since winter solstice']= 2
    d.loc[(d['Latitude'] <0) & (d['Month'] == 9),'months since winter solstice']= 3
    d.loc[(d['Latitude'] <0) & (d['Month'] == 10),'months since winter solstice']= 4
    d.loc[(d['Latitude'] <0) & (d['Month'] == 11),'months since winter solstice']= 5
    d.loc[(d['Latitude'] <0) & (d['Month'] == 12),'months since winter solstice']= 6
    return(d)



def volume_cone(d, h):
    volume = math.pi*((d/2)**2)*(h/3)
    return(volume)

def volume_sphere(d):
    volume = (1/6)*math.pi*(d**3)
    return(volume)

def morphometric_size_estimate(
        species,
        reference,
        coccosphere_length_min = None,
        coccosphere_length_max = None,
        coccosphere_width_min = None,
        coccosphere_width_max = None,
        coccolith_width_min = None,
        coccolith_width_max = None,
        coccolith_thickness = None,
        coccosphere = None,
        sample_n = 1000):

    coccolith_thickness_sd = None
    coccolith_thickness_mean = None
    coccosphere_length_sd = None
    coccosphere_width_sd = None
    cell_width = None

    if coccosphere_length_max is not None:
        coccosphere_length_sd = (coccosphere_length_max-coccosphere_length_min)/4
        coccosphere_length_mean = (coccosphere_length_max+coccosphere_length_min)/2
    
    if coccosphere_width_max is not None:
        coccosphere_width_sd = (coccosphere_width_max-coccosphere_width_min)/4
        coccosphere_width_mean = (coccosphere_width_max+coccosphere_width_min)/2

    if coccolith_width_max is not None:
        coccolith_thickness_sd = (coccolith_width_max-coccolith_width_min)/4
        coccolith_thickness_mean = (coccolith_width_max+coccolith_width_min)/2


    if (coccosphere_length_sd is not None) and (coccolith_thickness_sd is not None):
        cell_length = np.clip(np.random.normal(coccosphere_length_mean, coccosphere_length_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)

    if (coccosphere_width_sd is not None) and (coccolith_thickness_sd is not None):
        cell_width = np.clip(np.random.normal(coccosphere_width_mean, coccosphere_width_sd, sample_n), 0, None) - 2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None)



    if (coccolith_thickness is not None) and (coccosphere is None):
        cell_length = np.clip(np.random.normal(coccosphere_length_mean, coccosphere_length_sd, sample_n), 0, None) - 2*coccolith_thickness
    if (coccosphere_width_sd is not None) and (coccolith_thickness is not None) and (coccosphere is None):
        cell_width = np.clip(np.random.normal(coccosphere_width_mean, coccosphere_width_sd, sample_n), 0, None) - 2*coccolith_thickness



    if (coccolith_thickness is not None) and (coccosphere is not None):
        cell_length = coccosphere - 2*coccolith_thickness

    if (coccolith_thickness_sd is not None) and (coccosphere is not None):
        cell_length = coccosphere -(2*np.clip(np.random.normal(coccolith_thickness_mean, coccolith_thickness_sd, sample_n), 0, None))


    if cell_width is not None:
        cell_volume = volume_cone(cell_width, cell_length)
        shape = 'cone'
    else:
        cell_volume = volume_sphere(cell_length)
        shape = 'sphere'


    cell_mean = np.mean(cell_volume)
    try:
        cell_sd = np.sd(cell_volume)
    except:
        cell_sd = 0


    d = {'species':[species], 
        'mean': [cell_mean],
        'sd': [cell_sd],
        'shape':[shape],
        'ref': reference}

    d = pd.DataFrame(d)
    return(d)



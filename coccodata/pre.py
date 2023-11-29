
import pandas as pd
import numpy as np
import sys
import geopandas as gpd
from shapely.geometry import Point # Point class

from yaml import load, Loader
import glob, os


PATH_TO_SHAPEFILE = '/home/phyto/CoccoML/data/provinces/Longhurst_Biogeographical_Provinces.shp'

def preprocess_data(path):

    all_files = glob.glob(os.path.join(path, "*.csv"))

    d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)

    d = d.replace({0:pd.NA})

    d.reset_index(drop=False, inplace=True)


    references = pd.unique(d['Reference'])
    species = d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method']).columns

    for i in range(0, len(references)):

        mask = (d["Reference"] == references[i])
#        print(d[d["Reference"] == references[i]]["Method"])

        if d[d["Reference"] == references[i]]["Method"].iloc[0] == "SEM":
            if d[d["Reference"] == references[i]]["Reference"].iloc[0] != "Kleijne1984":
                d[species] = d[species].mask(mask, d[species].fillna(0))

        else:
            for j in range(0, len(species)):
                if d[d["Reference"] == references[i]][species[j]].sum()>0:
                    d[species[j]] = d[species[j]].mask(mask, d[species[j]].fillna(0))

                else:
                    None
    d['Month'] = d['Month'].astype('int64')
    d['Year'] = d['Year'].astype('int64')

    d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'], inplace=True)

    return(d)


d = preprocess_data("/home/phyto/CoccoData/data/abundances/")
d.rename(columns=lambda x: x.strip(), inplace=True)


with open('/home/phyto/CoccoData/data/classification/synonyms.yml', 'r') as f:
    groupings = load(f, Loader=Loader)

species = d.reset_index().drop(columns=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).columns

dict = {species:k
    for k, v in groupings.items()
    for species in v}

#df = d[species_observed]
df = d.copy()

df = (df.rename(columns=dict)
       .groupby(level=0, axis=1, dropna=False)).sum( min_count=1)

df.reset_index(inplace=True)

df = df.groupby(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).agg('mean')

df = df.drop(columns=['Phaeocystis pouchetii'])

counts =pd.DataFrame({'count': np.count_nonzero(df.fillna(0), axis=0), 'species': df.columns})
counts.to_csv("/home/phyto/CoccoData/raw_species_observations.csv" )

non_zero_spp = counts[counts['count']>0]['species']

df = df[non_zero_spp]
df.to_csv("/home/phyto/CoccoData/raw_observations.csv")

d = pd.read_csv("/home/phyto/CoccoData/raw_observations.csv")

depth_bins = np.linspace(-1, 300, 62).astype(np.int64) 
depth_labels = np.linspace(0, 300, 61).astype(np.int64) 
d = d[d["Depth"] >= 0]
d['Depth'] = pd.cut(d['Depth'], bins=depth_bins, labels=depth_labels).astype(np.int64) 

lat_bins = np.linspace(-90, 90, 181)
lat_labels = np.linspace(-89.5, 89.5, 180)
d['Latitude'] = pd.cut(d['Latitude'].astype(np.float64), bins=lat_bins, labels=lat_labels).astype(np.float64) 

lon_bins = np.linspace(-180, 180, 361)
lon_labels = np.linspace(-179.5, 179.5, 360)
d['Longitude'] = pd.cut(d['Longitude'].astype(np.float64), bins=lon_bins, labels=lon_labels).astype(np.float64) 

species = d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).columns

d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Month', 'Year'])[species].agg('mean').reset_index()
d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year'])

d.to_csv("/home/phyto/CoccoData/gridded_abundances.csv")

print("fin")

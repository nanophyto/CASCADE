
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


with open('/home/phyto/CoccoData/classification/synonyms.yml', 'r') as f:
    groupings = load(f, Loader=Loader)

species = d.reset_index().drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method']).columns

dict = {species:k
    for k, v in groupings.items()
    for species in v}

#df = d[species_observed]
df = d.copy()

df = (df.rename(columns=dict)
       .groupby(level=0, axis=1, dropna=False)).sum( min_count=1)

df.reset_index(inplace=True)


# df['Latitude'] = d['Latitude']
# df['Longitude'] = d['Longitude']
# df['Month'] = d['Month']
# df['Depth'] = d['Depth']
# df['Year'] = d['Year']
# df['Reference'] = d['Reference']
# df['Method'] = d['Method']


depth_bins = np.linspace(-1, 300, 62).astype(np.int64) 
depth_labels = np.linspace(0, 300, 61).astype(np.int64) 
df = df[df["Depth"] >= 0]

df['Depth'] = pd.cut(df['Depth'], bins=depth_bins, labels=depth_labels).astype(np.int64) 
df['Latitude'] = df['Latitude'].astype(np.float64).round()
df['Longitude'] = df['Longitude'].astype(np.float64).round()
df = df.groupby(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method']).agg('mean')


df = df.drop(columns=['Phaeocystis pouchetii'])

counts =pd.DataFrame({'count': np.count_nonzero(df.fillna(0), axis=0), 'species': df.columns})
counts.to_csv("/home/phyto/CoccoData/raw_species_observations.csv" )


non_zero_spp = counts[counts['count']>0]['species']


df = df[non_zero_spp]
df.to_csv("/home/phyto/CoccoData/raw_observations.csv")



d = pd.read_csv("/home/phyto/CoccoData/raw_observations.csv")


with open('/home/phyto/CoccoData/classification/synonyms.yml', 'r') as f:
    groupings = load(f, Loader=Loader)



#define depth bins    
bins = [-1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
        90, 95, 100, 125, 150, 175, 202, 225, 250, 275, 300]
labels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
          90, 95, 100, 125, 150, 175, 200, 225, 250, 275]

print(d.info(verbose=True))
d['Depth'] = pd.cut(d['Depth'], bins=bins, labels=labels).astype(np.int64) 


lat_bins = np.linspace(-90, 90, 181)
lat_labels = np.linspace(-89.5, 89.5, 180)
d['Latitude'] = pd.cut(d['Latitude'].astype(np.float64), bins=lat_bins, labels=lat_labels).astype(np.float64) 


lon_bins = np.linspace(-180, 180, 361)
lon_labels = np.linspace(-179.5, 179.5, 360)
d['Longitude'] = pd.cut(d['Longitude'].astype(np.float64), bins=lon_bins, labels=lon_labels).astype(np.float64) 



species = d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method']).columns

d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Month'])[species].agg('mean').reset_index()


dict = {species:k
    for k, v in groupings.items()
    for species in v}

df = (d.rename(columns=dict)
       .groupby(level=0, axis=1, dropna=False)).sum( min_count=1)



df['Latitude'] = d['Latitude']
df['Longitude'] = d['Longitude']
df['Month'] = d['Month']
df['Depth'] = d['Depth']

#df.set_index(['Latitude', 'Longitude', 'Depth', 'Month'], inplace=True)
df.to_csv("/home/phyto/CoccoData/abundances.csv", index=True)

d =df
d.set_index(['Latitude', 'Longitude', 'Depth', 'Month'], inplace=True)
d['dummy'] = 1


print("loading env")

df = pd.read_csv("/home/phyto/CoccoML/data/envdata_final.csv")
df = df[["temperature", "si", "phosphate", "din", "o2", "mld", "DIC", "irradiance", "chlor_a", "time", "depth", "lat", "lon"]]
df.rename(columns = {'depth':'Depth', 'lat':'Latitude', 'lon':'Longitude', 'time':'Month'}, inplace = True)

lat_bins = np.linspace(-90, 90, 181)
lat_labels = np.linspace(-89.5, 89.5, 180)
df['Latitude'] = pd.cut(df['Latitude'].astype(np.float64), bins=lat_bins, labels=lat_labels).astype(np.float64) 


lon_bins = np.linspace(-180, 180, 361)
lon_labels = np.linspace(-179.5, 179.5, 360)
df['Longitude'] = pd.cut(df['Longitude'].astype(np.float64), bins=lon_bins, labels=lon_labels).astype(np.float64) 


#define depth bins    
bins = [-1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
        90, 95, 100, 125, 150, 175, 202, 225, 250, 275, 300]
labels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
          90, 95, 100, 125, 150, 175, 200, 225, 250, 275]

#print(d.info(verbose=True))
df['Depth'] = pd.cut(df['Depth'], bins=bins, labels=labels).astype(np.int64) 



df = df.groupby(['Latitude', 'Longitude', 'Depth', 'Month']).agg('mean').reset_index()





df = df.set_index(['Latitude', 'Longitude', 'Depth', 'Month'])
df = df.dropna()


out = pd.concat([df, d], axis=1)

out = out[out["dummy"] == 1]
out = out.drop([ 'dummy'], axis = 1)


#define provinces:


out.reset_index(inplace=True)
#out.loc[out["Latitude"] >= 50, "Gephyrocapsa ericsonii"] = 0
#out.loc[out["Latitude"] <= -50, "Gephyrocapsa ericsonii"] = 0

counts =pd.DataFrame({'count': np.count_nonzero(out.fillna(0), axis=0), 'species': out.columns})
counts.to_csv("/home/phyto/CoccoData/final_species_observations.csv" )



# poly = gpd.read_file(PATH_TO_SHAPEFILE)

# for i in range(out.shape[0]):
#     for j in range(0, poly.shape[0]):
#         if poly['geometry'][j].contains(Point((out['Longitude'][i], out['Latitude'][i]))):
#             out.loc[i, 'FID'] = poly['ProvCode'][j]

# out = out.dropna(subset=['din'])


# out = out.set_index(['Latitude', 'Longitude', 'Depth', 'Month'])



# sys.path.insert(0, '/home/phyto/planktonSDM/functions/')
# from tune import tune 
# with open('/home/phyto/planktonSDM/configuration/devries2023_model_config.yml', 'r') as f:
#     model_config = load(f, Loader=Loader)
# model_config['remote'] = False

# traits = pd.read_csv("/home/phyto/CoccoML/data/traits.csv")
# species =  traits['species']


# metric = "shannon"
# measure = alpha_diversity(metric, out[species])
# out[metric] = measure.values


# metric = "simpson"
# measure = alpha_diversity(metric, out[species])
# out[metric] = measure.values

# metric = "observed_otus"
# measure = alpha_diversity(metric, out[species])
# out[metric] = measure.values



# out['total'] = out[species].sum( axis='columns')
# out['total_log'] = np.log(out['total'])


out.to_csv("/home/phyto/CoccoData/abundances_environment.csv", index=True)



# df = pd.read_csv("/home/phyto/CoccoML/data/envdata.csv")
# df = df[["temperature", "si", "phosphate", "din", "o2", "mld", "Fe", "DIC", "irradiance", "water_depth", "chlor_a", "time", "depth", "lat", "lon"]]


# df.rename(columns = {'depth':'Depth', 'lat':'Latitude', 'lon':'Longitude', 'time':'Month'}, inplace = True)


# lat_bins = np.linspace(-90, 90, 181)
# lat_labels = np.linspace(-89.5, 89.5, 180)
# df['Latitude'] = pd.cut(df['Latitude'].astype(np.float64), bins=lat_bins, labels=lat_labels).astype(np.float64) 


# lon_bins = np.linspace(-180, 180, 361)
# lon_labels = np.linspace(-179.5, 179.5, 360)
# df['Longitude'] = pd.cut(df['Longitude'].astype(np.float64), bins=lon_bins, labels=lon_labels).astype(np.float64) 



# df = df.set_index(['Latitude', 'Longitude', 'Depth', 'Month'])

# out = pd.concat([df, d], axis=1)

# out = out[out["dummy"] == 1]

# out = out[out["water_depth"] >= 200]
# out = out.drop(['water_depth', 'dummy'], axis = 1)

# print(pd.unique(out["FID"]))

# out.dropna(subset=['FID'], how='all', inplace=True)
# print(pd.unique(out["FID"]))

# out.to_csv("/home/phyto/CoccoML/data/abundances_environment.csv", index=True)


# counts =pd.DataFrame({'count': np.count_nonzero(out.fillna(0), axis=0), 'species': out.columns})
# counts.to_csv("/home/phyto/CoccoML/data/species_observations.csv" )


print("fin")


import pandas as pd
import numpy as np
import sys
import geopandas as gpd
from shapely.geometry import Point # Point class

from yaml import safe_load, load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from skbio.diversity import alpha_diversity

PATH_TO_SHAPEFILE = '/home/phyto/CoccoML/data/provinces/Longhurst_Biogeographical_Provinces.shp'

def preprocess_data():



    with open('/home/phyto/CoccoML/species_new.yml', 'r') as f:
        new_names = load(f, Loader=Loader)


    def clean_hagino200(): 
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Hagino2006.csv", na_values = "-")

        d['Latitude'] = d['Latitude'].str.rstrip('’S')
        d['Latitude'] = d['Latitude'].str.rstrip('’N')
        d['Latitude'] = d['Latitude'].str.replace('˚', '.')


        d['Longitude'] = d['Longitude'].str.rstrip('’E')
        d['Longitude'] = d['Longitude'].str.rstrip('’W')
        d['Longitude'] = d['Longitude'].str.replace('˚', '.')


        date = pd.to_datetime(d['Date'], format="%d/%m/%Y")

        d.rename(columns=new_names, inplace=True)

        d = d.groupby(by=d.columns, axis=1).sum()

        d['Reference'] = 'Hagino2006'
        d['Date'] = date
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns='Date', inplace=True)

        d = d.replace('--',np.NaN)

        d['Method'] = "SEM"
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d.fillna(0, inplace=True)

        return(d)


    def clean_devries2021(): 

        d = pd.read_csv("/home/phyto/CoccoData/abundances/deVries-etal_2020.csv")

        d['Date/Time'] = pd.to_datetime(d['Date/Time'])

        d = d[d.Reference != 'Takahashi & Okada (2000)'] #wrong data
        d = d[d.Reference != 'Godrijan et al. (2018)'] #shallow med
        d = d[d.Reference != 'Luan et al. (2016)'] #shallow china sea
        d = d[d.Reference != 'Cerino et al. (2017)'] #shallow med

        d.rename(columns = {'Depth water [m]':'Depth'}, inplace = True)

        #remove samples below 200m (9 samples)
        d = d[d.Depth <= 201] 


        #Saavedra-Pellitero et al. (2014) samples:
        #weirdly lablled lon

        d = d.replace(-182.73, 2.73)
        d = d.replace(-182.85, 2.85)
        d = d.replace(-185.47, 5.47)

        d.rename(columns=new_names, inplace=True)

        d = d.groupby(by=d.columns, axis=1).sum()

        d.loc[d["Reference"] == "Saavedra-Pellitero et al. (2014)", "Gephyrocapsa ericsonii"] = 0


        d['Month'] = pd.DatetimeIndex(d['Date/Time']).month
        d['Year'] = pd.DatetimeIndex(d['Date/Time']).year
        d.drop(columns=['Date/Time',"Sample method"], inplace=True)

        d.rename(columns = {'Depth water [m]':'Depth'}, inplace = True)

        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])


        return(d)


    def clean_okada1973(): 
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Okada_corrected.csv")
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d = d.round()
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()
    

        return(d)


    def clean_cortes2001():
        d = pd.read_csv("/home/phyto/CoccoData/abundances/cortes2001.csv", na_values = "ND")
        
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year

        d.drop(columns=['Date'], inplace=True)
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()
        

        return(d)

    def clean_takahashi2000():

        d = pd.read_csv("/home/phyto/CoccoData/abundances/Takahashi_raw.csv")

        counts =  d['counts']
        total = d['total']
        spp = d.drop(['Latitude', 'Longitude', 'Year', 'Depth', 'Month', 'counts', 'total'], axis = 1).reset_index(drop=True)


        for x in range(0, spp.shape[0]):
            for i in range(0, spp.shape[1]):
                spp.iloc[x][i] = np.round((spp.iloc[x][i]/counts[i]) * total[i])


        env = d[['Latitude', 'Longitude', 'Depth', 'Month', 'Year']]

        d =pd.concat([spp, env], axis=1) 
        d.fillna(0, inplace=True)
        d['Reference'] = "Takahashi2000"

        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])

        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()


        return(d)


    def clean_sal2013():

        d = pd.read_csv("/home/phyto/CoccoData/abundances/sal2013.csv")
        d['Reference'] = "Sal2013"
        d.rename(columns = {'Lat':'Latitude'}, inplace = True)
        d.rename(columns = {'Lon':'Longitude'}, inplace = True)
        d['Method'] = "LM"

        d = d[d['Latitude'].notna()]
        d = d[d['Longitude'].notna()]
        d = d[d['Depth'].notna()]
        d = d[d['Month'].notna()]
        d = d[d['Year'].notna()]


        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])


        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()



        return(d)

    def clean_obrien2013():


        with open('/home/phyto/CoccoData/references.yml', 'r') as f:
            references = load(f, Loader=Loader)

        #reverse key-binding to: long-short
        references = {v: k for k, v in references.items()}


        d = pd.read_csv("/home/phyto/CoccoData/abundances/obrien2013.csv")

        d = d.replace(["Utermohl", 
                "phase-contrast microscope",
                "Sedimentation and Inverted Microscope",
                    "Inverted Microscope",
                    "Light Microscope (x400)",
                    "epifluorescence microscopy",
                    "inverted microscope",
                    "light microscopy",
                "Inverted light and epifluorescent microscopy (Utermohl)",
                'Inverted microscope (Utermohl)',
                'Inverted microscope',
                ], "LM")

        print(d['Instrument/Method'].unique())

        d = d[d["Instrument/Method"].isin(["LM"])]

        d.rename(columns = {'Method References':'Reference'}, inplace = True)
        d.rename(columns = {'Depth (m)':'Depth'}, inplace = True)
        d.rename(columns = {'Longitutude':'Longitude'}, inplace = True)
        d.rename(columns = {'Instrument/Method':'Method'}, inplace = True)


        d = d.pivot_table(index=['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'], columns=['Species/lowest classification'], values=['Cells/l'])
        d.columns = d.columns.droplevel(0)
        d = d.reset_index().rename_axis(None, axis=1)

        #d.reset_index(inplace=True)
        d=d.replace({"Reference": references})
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])

        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()
        
        return(d)



    def clean_estrada2016():
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Malaspina2010.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year

        d.drop(columns=['Date', 'Station'], inplace=True)

        d['Method'] = "LM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()
        
        return(d)


    def clean_baumann2000():
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Baumann2000.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year

        d.drop(columns=['Date'], inplace=True)
        d['Reference'] = "Baumann2000"
        d['Method'] = "SEM"


        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()

        return(d)

    def clean_kleijne1984():
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Kleijne1984.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns=['Date'], inplace=True)
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()

        return(d)

    def clean_keuter2023():
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Keuter2023.csv")
        d['Method'] = "SEM"
        d['Reference'] = "Keuter2023"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()

        return(d)
    
    def clean_keuter2022():
        d = pd.read_csv("/home/phyto/CoccoData/abundances/Keuter2022.csv")
        d['Method'] = "SEM"
        d['Reference'] = "Keuter2022"

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns=['Date'], inplace=True)

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000
        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()

        return(d)
    

    estrada2016 = clean_estrada2016()
    hagino2006 = clean_hagino200()
    devries2021 = clean_devries2021()
    okada1973 = clean_okada1973()
    takahashi2000 = clean_takahashi2000()

    cortes2001 = clean_cortes2001()
    baumann2000 = clean_baumann2000()
    sal2013 = clean_sal2013()
    obrien2013 = clean_obrien2013()
    kleijne1984 = clean_kleijne1984()
    keuter2023 = clean_keuter2023()
    keuter2022 = clean_keuter2022()

    d = pd.concat([hagino2006, devries2021, okada1973, takahashi2000, estrada2016, cortes2001, 
                   baumann2000, sal2013, obrien2013, kleijne1984, keuter2023, keuter2022])
#    print(d.head())

#    with np.printoptions(threshold=np.inf):
#        print(d.columns)


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


d = preprocess_data()
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

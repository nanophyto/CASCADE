import pandas as pd
import numpy as np
from yaml import  load, Loader
from pandas.api.types import is_integer_dtype
from datetime import datetime

#load raw observations:
d = pd.read_csv("/home/phyto/CoccoData/raw_observations.csv")

"""

Abundance dataset quality control:

Coordinates make sense

Species are defined

Species are non-zero 

References are in library (TO DO)

Mean values match those reported in literature (TO DO)

Calculations and generation of table for paper:
    - number of observations
    - time period covered
    - mean abundance
    - reference
    - provinces included in data set (optional - if there is time)

"""


class ValueCheck():
    """

    Notes
    -----
    Assumes columns are named: 
    'Depth', 'Latitude', 'Longitude', 'Month', 'Year', 'Reference', 'Method'

    """
    def __init__(self, d):
        self.d = d

    def coords_check(self):

        #check if longitude falls within range -180 to 180
        if (self.d['Longitude'].any()<-180) or (self.d['Longitude'].any()>180):
            raise ValueError("Longitude outside of range -180 to 180")

        if not (
            (self.d['longitude'] > -180).all() &
            (self.d['longitude'] < 180).all()
        ):
            raise ValueError("Longitude outside...")


        #check if latitude falls within range -90 to 90
        if (self.d['Latitude'].any()<-90) or (self.d['Latitude'].any()>90):
            raise ValueError("Latitude outside of range -90 to 90")

        #check if depth is postitive
        if (self.d['Depth'].any()<0):
            raise ValueError("Depth not stricly positive")

        #check if month is an interger and falls within range 1 to 12
        if self.d['Month'].dtypes != 'int64':
            raise ValueError("Month is not an integer")
        if not 1 <= self.d['Month'].any() <= 12:
            raise ValueError("Month is greater than 12 or less than 1")

        #check if year is an interger
        if self.d['Year'].dtypes != 'int64':
            raise ValueError("Year is not an integer")
        
        #check if year is in: 1900-now
        now = datetime.now().year
        if 1900 > self.d['Year'].any() > now:
            raise ValueError("Year does not fall between: 1900-" + str(now))

    def methods_check(self, methods=None):
        if methods != None:
            #check if method is defined:
            obs = self.d['Method'].unique()
            if not obs.any() in methods:
                raise ValueError("data contains undefined method")
        else:
            raise ValueError("data contains undefined method")



    #check if reference is in reference library
    def reference_check(self, methods=None):
        with open('/home/phyto/CoccoML/refences.yml', 'r') as f:
            groupings = load(f, Loader=Loader)



    def species_check(self, species_yaml):

        #set lat, lon, depth, month, year, method as index
        self.d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'], inplace=True)

        #rename and group synonyms 
        with open(species_yaml, 'r') as f:
            groupings = load(f, Loader=Loader)

        self.d.rename(columns=lambda x: x.strip(), inplace=True)

        #drop index
        species = self.d.reset_index().drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method']).columns
        #could be more efficient

        dict = {species:k
            for k, v in groupings.items()
            for species in v['alt']}

        df = self.d[species]

        df = (df.rename(columns=dict)
            .groupby(level=0, axis=1, dropna=False)).sum( min_count=1)

        #check if final species are in species library
        species_library = {v: k for k, v in dict.items()}

        species_observed = df.columns

        if (set(species_observed).issubset(species_library)):
            print("all species are defined")
        else:
            raise ValueError("undefined species:" + str(set(species_observed).difference(species_library)))

        #check if abundances are positive
        if not (self.d.min(axis=1).any() > 0).all():
            raise ValueError("negative abundance values")



    #print number of observations:
#    print("number of observations: "+ str(np.sum(df.count())))


    def regrid_data(self):

        species = self.d.reset_index().drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method']).columns
        depth_bins = np.linspace(0, 205, 62)
        depth_labels = np.linspace(0, 200, 61)
        self.d['Depth'] = pd.cut(self.d['Depth'], bins=depth_bins, labels=depth_labels).astype(np.int64) 

        lat_bins = np.linspace(-90, 90, 181)
        lat_labels = np.linspace(-89.5, 89.5, 180)
        self.d['Latitude'] = pd.cut(self.d['Latitude'].astype(np.float64), bins=lat_bins, labels=lat_labels).astype(np.float64) 

        lon_bins = np.linspace(-180, 180, 361)
        lon_labels = np.linspace(-179.5, 179.5, 360)
        self.d['Longitude'] = pd.cut(self.d['Longitude'].astype(np.float64), bins=lon_bins, labels=lon_labels).astype(np.float64) 


        self.d_sd = self.d.groupby(['Latitude', 'Longitude', 'Depth', 'Month'])[species].agg('sd').reset_index()
        self.d = self.d.groupby(['Latitude', 'Longitude', 'Depth', 'Month'])[species].agg('mean').reset_index()


m = ValueCheck(d)
m.coords_check()
m.methods_check(methods=['SEM','LM'])
m.species_check(species_yaml='/home/phyto/CoccoData/groupings.yml')















d = pd.read_csv("/home/phyto/CoccoML/data/raw_observations.csv")
print(len(d['Reference'].unique()))


#check if all genera are defined in sheward data set:
genera = pd.read_csv("/home/phyto/CoccoML/data/species_genera.csv")

sheward = pd.read_csv("//home/phyto/CoccoML/data/PICPOC_raw_V2.csv")

sheward['Genus'] = sheward['Genus'].str.strip()

if (set(genera['genera']).issubset(sheward['Genus'])):
    print("all genera are defined")
else:
    raise ValueError("undefined genera:" + str(set(genera['genera']).difference(sheward['Genus'])))







# out = out.dropna(subset=['din'])




#estimate size


#estimate POC


#estimate PIC


#monthly climatology


#export



#appendix

print("fin")
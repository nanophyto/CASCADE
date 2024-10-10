import pandas as pd
import numpy as np
from yaml import load, Loader
import glob, os
from functions import rename_synonyms, rename_synonyms_wide, morphometric_size_estimate, volume_sphere
import math 

class pre_allometry():

    def __init__(self, import_path = "../data/unprocessed/allometry/", 
                 export_path = "../data/allometry/",
                 classification_path = "../data/classification/"):
        self.import_path = import_path
        print(self.import_path)
        self.export_path = export_path
        self.classification_path = classification_path

    def pre_sheward_allometry(self):
        d = pd.read_csv(self.import_path + "sheward2024.csv")
        d = d.dropna()
        d = d[d['PIC pg C'] >0]
        d = d[d['Volume'] >0]
        d['species'] = d['Genus'] + " " + d['Species']
        d = d.rename(columns={"Volume": "volume", "PIC pg C": "pg pic", "Phase":"phase"})
        d = d[['species', 'volume', 'pg pic', 'phase']]
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml', remove_duplicate=False)
        filter = d['species'].str.contains('undefined')
        d = d[~filter]
        d.to_csv(self.export_path + "sheward2024.csv", index=False)
        print("finished processing: sheward2024")


class pre_size():

    def __init__(self, import_path = "../data/unprocessed/sizes/", 
                 export_path = "../data/sizes/",
                 classification_path = "../data/classification/"):
        self.import_path = import_path
        self.export_path = export_path
        self.classification_path = classification_path


    def resample_size(self, d, species):
        d = d[d["species"]==species]


        if len(d)>1:
            size = []
            for i in range(len(d)):
                size_distribution = np.random.normal(d.iloc[i]['mean'], d.iloc[i]['sd'], 10000)
                size.append(size_distribution)
            mean = np.nanmean(size)
            sd = np.nanstd(size)
        else:
            mean = d['mean']
            sd = d['sd']

        species_data = {'species': [species],
                    'mean': mean,
                    'sd': sd}
        
        species_data = pd.DataFrame(species_data)

        return(species_data)

    def resample_size_d(self, d):
        species_list = d['species'].unique()

        new_d = []

        for i in range(len(species_list)):
            size = self.resample_size(d, species_list[i])
            new_d.append(size)

        d = pd.concat(new_d)

        return(d)

    def pre_villiot2021a(self):    
        d = pd.read_csv(self.import_path + "viliot2021a.csv")
        d = d.rename(columns={'std': "sd"})
        d = self.resample_size_d(d)
        d = d[['species', 'mean', 'sd']] 
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d['reference'] = "villiot2021a"
        d['method'] = 'light microscopy'
        d['n'] = 100
        d.to_csv(self.export_path + "viliot2021a.csv", index=False)

    def pre_villiot2021b(self):
        d = pd.read_csv(self.import_path + "villiot2021b.csv")
        d = d.rename(columns={'cell volume': "mean"})
        d['mean'] = np.round(d['mean'])
        d['sd'] = None
        d = d[['species', 'mean', 'sd']] 
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d['reference'] = "villiot2021b"
        d['method'] = 'literature morphometrics'
        d['n'] = 1
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d.to_csv(self.export_path + "viliot2021b.csv", index=False)

    def pre_obrien2013a(self):
        d = pd.read_csv(self.import_path + "obrien2013.csv")
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d = d.rename(columns={'diameter': "mean"})
        d['mean'] = (1/6)*math.pi*(d['mean']**3)
        d['mean'] = np.round(d['mean'])
        d['reference'] = 'obrien2013a'
        d['sd'] = None
        d['method'] = 'light microscopy'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "obrien2013a.csv", index=False)
        return(d)

    def pre_devries2024(self):

        aurisinae = morphometric_size_estimate(
                species = 'Syracosphaera aurisinae',
                reference = "lecal1951",
                coccosphere_length_min = 23,
                coccosphere_length_max = 26,
                coccosphere_width_min = 13,
                coccosphere_width_max = 17,
                coccolith_width_min = 1.1,
                coccolith_width_max = 1.2)

        caudatus = morphometric_size_estimate(
                species = 'Calciopappus caudatus',
                reference = "gaarder1954",
                coccosphere_length_min = 26,
                coccosphere_length_max = 36,
                coccosphere_width_min = 3.5,
                coccosphere_width_max = 4,
                coccolith_thickness = 1)

        blokii = morphometric_size_estimate(
            species = 'Calicasphaera blokii',
            reference = 'cros2002',
            coccosphere = 6,
            coccolith_thickness = 1)

        concava = morphometric_size_estimate(
            species = 'Calicasphaera concava',
            reference = 'cros2002',
            coccosphere = 6,
            coccolith_thickness = 1.3)

        borealis = morphometric_size_estimate(
                species = 'Syracosphaera borealis',
                reference = 'okada1977',
                coccosphere_length_min = 6.5, 
                coccosphere_length_max = 8.2,
                coccolith_thickness = 1)


        adenensis = morphometric_size_estimate(
                species = 'Sphaerocalyptra adenensis',
                reference = 'cros2002',
                coccosphere_length_min = 5.5, 
                coccosphere_length_max = 8.5,
                coccolith_thickness = 1)


        gaudii_POL = morphometric_size_estimate(
                species = 'Alisphaera gaudii POL',
                reference = 'cros2002',
                coccosphere_length_min = 5.6, 
                coccosphere_length_max = 10.6,
                coccolith_thickness = 1)


        pienaarii = morphometric_size_estimate(
                species = 'Helladosphaera pienaarii',
                reference = ['norris1985, kleijne1991'],
                coccosphere_length_min = 6, 
                coccosphere_length_max = 11,
                coccolith_thickness = 0.9)

        strigilis = morphometric_size_estimate(
                species = 'Syracosphaera strigilis',
                reference = 'cros2002',
                coccosphere_length_min = 6, 
                coccosphere_length_max = 9,
                coccolith_thickness = 1)

        formosus = morphometric_size_estimate(
                species = 'Ophiaster formosus',
                reference = 'cros2002',
                coccosphere_length_min = 4.5, 
                coccosphere_length_max = 7.5,
                coccolith_width_min = 0.7,
                coccolith_width_max = 0.8)

        type_5 = morphometric_size_estimate(
                species = 'Pappomonas sp. type 5',
                reference = 'cros2002',
                coccosphere_length_min = 6,
                coccosphere_length_max = 7,
                coccolith_thickness = 1)

        squamosa = morphometric_size_estimate(
                species = 'Syracosphaera squamosa',
                reference = 'kleijne2009',
                coccosphere = 6,
                coccolith_width_min = 0.9,
                coccolith_width_max = 1.3)

        reniformis = morphometric_size_estimate(
                species = 'Syracosphaera reniformis',
                reference = 'kleijne2009',
                coccosphere_length_min = 6, 
                coccosphere_length_max = 8,
                coccolith_width_min = 1,
                coccolith_width_max = 1.2)

        sphaeroidea_hol = morphometric_size_estimate(
                species = 'Calyptrosphaera sphaeroidea HOL',
                reference = 'cros2002',
                coccosphere_length_min = 5.5, 
                coccosphere_length_max = 12,
                coccolith_thickness = 1)

        cristatus = morphometric_size_estimate(
                species = 'Ceratolithus cristatus',
                reference = ['obrien2013 and young1998'],
                coccosphere_length_min =  7,
                coccosphere_length_max = 18.9,
                coccolith_width_min = 0.2, 
                coccolith_width_max = 0.4 )

        confusus = morphometric_size_estimate(
                species = 'Helicosphaera HOL confusus type',
                reference = 'cros2002',
                coccosphere_length_min = 9,  
                coccosphere_length_max = 14,
                coccolith_thickness = 1 )

        arctica =  pd.DataFrame({'species':['Wigwamma antarctica'], 
                'mean': [volume_sphere(5)],
                'sd': [None],
                'shape':['sphere'],
                'ref': ['thomsen1988']})

        marsilii = morphometric_size_estimate(
                species = 'Zygosphaera marsilii',
                reference = 'cros2002',
                coccosphere_length_min = 6.5, 
                coccosphere_length_max = 8.5, 
                coccolith_thickness = 0.3)

        hydroideus = morphometric_size_estimate(
                species = 'Ophiaster hydroideus',
                reference = 'cros2002',
                coccolith_width_min = 0.7, 
                coccolith_width_max = 0.9, 
                coccosphere = 6 )

        d = pd.concat([type_5, aurisinae, sphaeroidea_hol,
            reniformis, squamosa, arctica, pienaarii, cristatus, hydroideus,
            marsilii, caudatus, borealis, adenensis, blokii, concava, gaudii_POL,
            strigilis, formosus, confusus])

        d.to_csv(self.import_path + "devries2024.csv", index=False) #export to unprocessed

        d = d[['species', 'mean', 'sd']] 
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d['reference'] = 'devries2024'
        d['method'] = 'literature morphometrics'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 1
        d.to_csv(self.export_path + "devries2024.csv", index=False)

    def pre_sheward2024(self):

        #read sheward size data
        d = pd.read_csv(self.import_path + "sheward2024.csv")

        d = d[d['Spheres Y/N?']!="Flattened"]
        d['Species'] = d['Species'].str.strip()
        d['Genus'] = d['Genus'].str.strip()
        d['species'] = d['Genus'] + " " + d['Species'] 

        d = d[['species', 'Estimated cell volume', 'PIC pg C']]

        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml', 
                            remove_duplicate=False, check_synonyms=False)

        d = d.groupby(by="species").agg(["mean", "std", "count"]).reset_index()

        d['sd'] = np.round(d['Estimated cell volume']['std'], 1)
        d['mean'] = np.round(d['Estimated cell volume']['mean'], 1)
        d['n'] = d['Estimated cell volume']['count']

        d = d[['species', 'sd', 'mean', 'n']]

        #rename because column names are tuples for some reason??
        d.columns = ['species', 'sd', 'mean', 'n']

        d['reference'] = 'sheward2024'
        d['method'] = "AMT morphometrics"
        d = d[['species', 'mean', 'sd', 'method', 'n', 'reference']] 
        d = d[d['n']>0]
        d.to_csv(self.export_path + "sheward2024.csv", index=False)


    def pre_young2024(self):
        d = pd.read_csv(self.import_path + "young2024.csv")
        d['mean'] = d['avg(volume)']
        d['sd'] = d['stddev(volume)']
        d['n'] = d['count(*)']
        d = d[['species', 'mean', 'sd', 'n']] 
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d['reference'] = 'young2024'
        d['method'] = 'light microscopy'
        d.to_csv(self.export_path + "young2024.csv", index=False)

    def pre_gafar2019(self):
        d = pd.read_csv(self.import_path + "gafar2019.csv")
        d['mean'] = d['cell volume']
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d = self.resample_size_d(d)
        d['reference'] = 'gafar2019'
        d['method'] = 'light microscopy'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "gafar2019.csv", index=False)

    def pre_fiorini2011(self):
        d = pd.read_csv(self.import_path + "fiorini2011.csv")
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d = self.resample_size_d(d)
        d['reference'] = 'fiorini2011'
        d['method'] = 'coulter counter'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "fiorini2011.csv", index=False)

    def pre_gerecht2018(self):
        d = pd.read_csv(self.import_path + "gerecht2018.csv")
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d = self.resample_size_d(d)
        d['reference'] = 'gerecht2018'
        d['method'] = 'light microscopy'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "gerecht2018.csv", index=False)

    def pre_oviedo2014(self):
        d = pd.read_csv(self.import_path + "oviedo2014.csv")
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d = self.resample_size_d(d)
        d['reference'] = 'oviedo2014'
        d['method'] = 'coulter counter'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "oviedo2014.csv", index=False)

    def pre_supraha2015(self):
        d = pd.read_csv(self.import_path + "supraha2015.csv")
        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml')
        d = self.resample_size_d(d)
        d['reference'] = 'supraha2015'
        d['method'] = 'light microscopy'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "supraha2015.csv", index=False)

    def preprocess_all(self):
        self.pre_villiot2021a()
        self.pre_villiot2021b()
        self.pre_obrien2013a()
        self.pre_devries2024()
        self.pre_sheward2024()
        self.pre_young2024()
        self.pre_gafar2019()
        self.pre_fiorini2011()
        self.pre_gerecht2018()
        self.pre_oviedo2014()
        self.pre_supraha2015()

        print("finished processing:")
        print("villiot2021a")
        print("villiot2021b")
        print("obrien2013a")
        print("devries2024")
        print("sheward2024")
        print("young2024")
        print("gafar2019")
        print("fiorini2011")
        print("gerecht2018")
        print("oviedo2014")
        print("supraha2015")

class pre_pic():

    def __init__(self, import_path = "../data/unprocessed/pic/", 
                 export_path = "../data/pic/",
                 classification_path= "../data/classification/"):
        self.import_path = import_path
        self.export_path = export_path
        self.classification_path = classification_path

    def pre_devries2024(self):
        def estimate_pontosphaera():
            pic_per_lith = [77.97, 133.07, 237.68] #Guerreiro2021
            observations_per_estimate = [13, 5, 2] #Guerreiro2021
            number_of_liths = [45]   #yang and wei 2003

            p=observations_per_estimate/np.sum(observations_per_estimate) 

            pic_sphere = np.random.choice(pic_per_lith, 10000, p=p)*np.random.choice(number_of_liths, 10000)

            mean_pic_per_sphere =  np.round(np.mean(pic_sphere), -2)*0.12
            sd_pic_per_sphere = np.round(np.std(pic_sphere), -2)*0.12

            d = {'species':['undefined Pontosphaera'], 
                'mean': [mean_pic_per_sphere],
                'sd': [sd_pic_per_sphere],
                'ref': ['guerreiro2021 and yang2003']}

            d = pd.DataFrame(d)

            return(d)


        def estimate_flabellatus():

            pic_per_lith = [9.65, 11.62, 4.25] #Guerreiro2021
            observations_per_estimate = [213, 106, 108] #Guerreiro2021
            number_of_liths = [40, 86, 56, 18, 120]   #yang and wei 2003

            p=observations_per_estimate/np.sum(observations_per_estimate) 

            pic_sphere = np.random.choice(pic_per_lith, 10000, p=p)*np.random.choice(number_of_liths, 10000)

            mean_pic_per_sphere =  np.round(np.mean(pic_sphere), -1)*0.12
            sd_pic_per_sphere = np.round(np.std(pic_sphere), -1)*0.12


            d = {'species':['Gladiolithus flabellatus'], 
                'mean': [mean_pic_per_sphere],
                'sd': [sd_pic_per_sphere],
                'ref': ['guerreiro2021 and yang2003']}

            d = pd.DataFrame(d)

            return(d)

        d_flabellatus = estimate_flabellatus()
        d_pontosphaera = estimate_pontosphaera()

        d = pd.concat([d_flabellatus, d_pontosphaera])
        d['reference'] = 'devries2024'
        d['method'] = "literature morphometrics"
        d['n'] = 1
        d.to_csv(self.export_path + "devries2024.csv", index=False)


    def pre_sheward2014(self):

        def cell_pic(cn, cl, ks):
            pic = (2.7*ks*cn*(cl**3))*0.12
            return(pic)

        df = pd.read_csv(self.import_path + "sheward2014.csv")
        df.replace({'Coccolithus braarudii': 'Coccolithus pelagicus'}, inplace=True)   

        def estimate_pic(d, species, ks):
            d = d[d['species']==species]
            d = df.assign(pic=cell_pic(d['CN'], d['CL'], ks)) 
            return(d[['pic', 'species']])    

        d = estimate_pic(df, 'Coccolithus pelagicus', ks=0.06)

        d = d.groupby(by="species").agg(["mean", "std", "count"]).reset_index()

        d['sd'] = np.round(d['pic']['std'], 1)
        d['mean'] = np.round(d['pic']['mean'], 1)
        d['n'] = d['pic']['count']

        d = d[['species', 'sd', 'mean', 'n']]
        d.reset_index(inplace=True, drop=True)
        d = d.dropna()

        d['reference'] = 'sheward2014'
        d['method'] = "lab morphometrics"

        df = pd.DataFrame({
            'reference' : d['reference'], 
            'method' : d['method'],
            'species': d['species'],
            'mean': d['mean'],
            'sd': d['sd'],
            'n': d['n'],
        })
        df.to_csv(self.export_path + "sheward2014.csv", index=False)


    def pre_sheward2016(self):

        def cell_pic(cn, cl, ks):
            pic = (2.7*ks*cn*(cl**3))*0.12
            return(pic)

        def estimate_leptoporus():
            d2 = pd.read_csv(self.import_path + "sheward2017_leptoporus.csv")
            d2['pic'] = cell_pic(d2['CN'], d2['CL'], 0.08)
            d2.drop(columns=['species'], inplace=True)

            d1 = pd.read_csv(self.import_path + "sheward2017_quadriperforatus.csv")
            d1['pic'] = cell_pic(d1['CN'], d1['CL'], 0.08)
            d1.drop(columns=['species'], inplace=True)

            d = pd.concat([d1, d2])
            d['species'] = 'Calcidiscus leptoporus'
            return(d[['pic', 'species']])

        def estimate_carteri():
            d = pd.read_csv(self.import_path + "sheward2017_carteri.csv")
            d['pic'] = cell_pic(d['CN'], d['CL'], 0.05)
            return(d[['pic', 'species']])

        d = pd.concat([estimate_leptoporus(), estimate_carteri()])

        d = d.groupby(by="species").agg(["mean", "std", "count"]).reset_index()

        d['sd'] = np.round(d['pic']['std'], 1)
        d['mean'] = np.round(d['pic']['mean'], 1)
        d['n'] = d['pic']['count']
        d = d[['species', 'sd', 'mean', 'n']]

        d['reference'] = 'sheward2016'
        d['method'] = "lab morphometrics"
        d.reset_index(inplace=True, drop=True)
        d = d.dropna()

        df = pd.DataFrame({
            'reference' : d['reference'], 
            'method' : d['method'],
            'species': d['species'],
            'mean': d['mean'],
            'sd': d['sd'],
            'n': d['n'],
        })

        df.to_csv(self.export_path + "sheward2016.csv", index=False)


    def pre_sheward2024(self):

        #read sheward size data
        d = pd.read_csv(self.import_path + "sheward2024.csv")

        d = d[d['Coccosphere condition']=="Intact"]
        d['Species'] = d['Species'].str.strip()
        d['Genus'] = d['Genus'].str.strip()
        d['species'] = d['Genus'] + " " + d['Species'] 

        d = d[['species', 'pg pic']]

        d = rename_synonyms(d, classification_path  = self.classification_path + 'synonyms.yml', 
                            remove_duplicate=False)

        d = d.groupby(by="species").agg(["mean", "std", "count"]).reset_index()
        d['sd'] = np.round(d['pg pic']['std'], 1)
        d['mean'] = np.round(d['pg pic']['mean'], 1)
        d['n'] = d['pg pic']['count']

        d = d[['species', 'sd', 'mean', 'n']]

        #rename because column names are tuples for some reason??
        d.columns = ['species', 'sd', 'mean', 'n']

        d['reference'] = 'sheward2024'
        d['method'] = "AMT morphometrics"
        d = d[['species', 'mean', 'sd', 'n', 'method', 'reference']] 
        d.to_csv(self.export_path + "sheward2024.csv", index=False)


    def preprocess_all(self):
        self.pre_sheward2024()
        self.pre_sheward2016()
        self.pre_sheward2014()
    #    self.pre_devries2024()

        print("finished processing:")
        print("sheward2014")
        print("sheward2016")
        print("sheward2024")
        print("devries2024")


class pre_abundances():

    def __init__(self, root,
                    import_path = "../data/unprocessed/allometry/", 
                    export_path = "../data/allometry/",
                    reference_path = "../data/references/",
                    classification_path = '../data/classification/synonyms.yml'):
        
        self.root = root
        print(root)
        self.import_path = import_path
        self.export_path = export_path
        self.yaml_path = reference_path
        self.classification_path = classification_path

    def clean_hagino200(self): 
        d = pd.read_csv(self.import_path + "Hagino2006.csv", 
                        na_values = "-")

        d['Latitude'] = d['Latitude'].str.rstrip('’S')
        d['Latitude'] = d['Latitude'].str.rstrip('’N')
        d['Latitude'] = d['Latitude'].str.replace('˚', '.')


        d['Longitude'] = d['Longitude'].str.rstrip('’E')
        d['Longitude'] = d['Longitude'].str.rstrip('’W')
        d['Longitude'] = d['Longitude'].str.replace('˚', '.')

        date = pd.to_datetime(d['Date'], format="%d/%m/%Y")

        #d = d.groupby(by=d.columns, axis=1).sum()
        d = d.T.groupby(level=0).sum().T

        d.insert(0, 'Reference', 'Hagino2006')

        d['Date'] = date
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns='Date', inplace=True)

        d = d.replace('--',np.NaN)
        d.insert(0, 'Method', 'SEM')

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.fillna(0, inplace=True)

        d = d.reset_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        
        d = rename_synonyms_wide(d, classification_path=self.classification_path)
        d.to_csv(self.export_path + "Hagino2006.csv", index=False)



    def clean_devries2021(self): 
        d = pd.read_csv(self.import_path + "devries2020.csv")

        d['Date/Time'] = pd.to_datetime(d['Date/Time'])

        d = d[d.Reference != 'Takahashi & Okada (2000)'] #wrong data
        d = d[d.Reference != 'Dimiza et al. (2015)'] #double entry in database

        d.rename(columns = {'Depth water [m]':'Depth'}, inplace = True)

        d.columns = d.columns.str.replace(' [#/l]', '')

        d = d.replace(-182.73, 2.73)
        d = d.replace(-182.85, 2.85)
        d = d.replace(-185.47, 5.47)

        d['Day'] = pd.DatetimeIndex(d['Date/Time']).day
        d['Month'] = pd.DatetimeIndex(d['Date/Time']).month
        d['Year'] = pd.DatetimeIndex(d['Date/Time']).year

        d.drop(columns=['Date/Time',"Sample method"], inplace=True)

        d.rename(columns = {'Depth water [m]':'Depth'}, inplace = True)

        d.insert(0, 'Method', 'SEM')

        with open(self.root + '/data/references/refs_devries.yml', 'r') as f:
            short_refs = load(f, Loader=Loader)

        d=d.replace({"Reference": short_refs})

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.reset_index(inplace=True)

        for i in range(len(d['Reference'].unique())):
            reference_to_export = d['Reference'].unique()[i]
            d[d.Reference==reference_to_export].to_csv(self.export_path + reference_to_export + ".csv", index=False)



    def clean_johnson2023(self):
        d = pd.read_csv(self.import_path + "johnson2023.csv")

        with open(self.yaml_path + 'johnson2023.yml', 'r') as f:
            references = load(f, Loader=Loader)

        d = d[d.Reference != 'Dimiza et al. 2016'] #double entry in database

        d=d.replace({"Reference": references})

        d = d.replace(["Scanning electron microscope (SEM)",
                    "SEM and LM",
                    "Polarized light microscopy (PLM) + Scanning electron microscope (SEM)",   
                ], "SEM")
        
        d = d.replace(["Light microscopy" 
                ], "LM")
        d.columns = [col.replace(' [#/l]', '') for col in d.columns]        
        d.drop(columns=['Date/Time'], inplace=True)
        d.rename(columns = {'Sample method':'Method'}, inplace = True)
        d.columns = [col.strip() for col in d.columns]

        d = d[d['Reference']!="Estrada1985"] #unpublished data (did not reach out for permission)

#        d = rename_synonyms_wide(d, classification_path=self.classification_path)

        for i in range(len(d['Reference'].unique())):
            reference_to_export = d['Reference'].unique()[i]
            d[d.Reference==reference_to_export].to_csv(self.export_path + reference_to_export + ".csv", index=False)



    def clean_okada1973(self): 

        d1 = pd.read_csv(self.import_path + "okada/Okada1973.csv")
        d2 = pd.read_csv(self.import_path + "okada/Okada1973B.csv")
        d3 = pd.read_csv(self.import_path + "okada/traverse2.csv")
        d4 = pd.read_csv(self.import_path + "okada/traverse4.csv")
        d5 = pd.read_csv(self.import_path + "okada/traverse4_A.csv")
        d6 = pd.read_csv(self.import_path + "okada/traverse5.csv")
        d7 = pd.read_csv(self.import_path + "okada/traverse13.csv")

        d = pd.concat([d1, d2, d3, d4, d5, d6, d7])

        #d.insert(0, 'Method', 'SEM')
        #d.insert(0, 'Reference', 'Okada1973')
        d = d.assign(Reference='Okada1973')
        d = d.assign(Method='SEM')


        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.drop(['Temperature', 'station', 'Sample', 'Date'], axis = 1)

        def cell_per_litre(d):

            #find the number of cells counted:
            total = d['Total'].values

            d = d.drop(['Total'], axis = 1)
            counts= d.sum(axis=1).values

            #loop over each row and estimate the cells per liter by dividing by the number of counts and multiplying by the total:
            
            for x in range(0, d.shape[0]): #for each row:
                for i in range(0, d.shape[1]): #for each column:
                    d.iloc[x, i] = np.round((d.iloc[x, i]/counts[x]) * total[x])
            return(d)

        d = cell_per_litre(d)
        d = d.reset_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.to_csv(self.export_path + "Okada1973.csv", index=False)



    def clean_cortes2001(self):
        d = pd.read_csv(self.import_path + "cortes2001.csv", na_values = "ND")
        
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date'], inplace=True)
        d.insert(0, 'Method', 'SEM')

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv(self.export_path + "Cortes2001.csv", index=False)


    def clean_takahashi2000(self):

        d = pd.read_csv(self.import_path + "takahashi2010.csv")
        d.insert(0, 'Reference', "Takahashi2000")
        d.insert(0, 'Method', 'SEM')
        d['Date/Time'] = pd.to_datetime(d['Date/Time'])
        d['Day'] = pd.DatetimeIndex(d['Date/Time']).day
        d['Month'] = pd.DatetimeIndex(d['Date/Time']).month
        d['Year'] = pd.DatetimeIndex(d['Date/Time']).year

        counts =  d['counts'].values
        total = d['total'].values

        d = d.drop(['counts', 'total', 'Date/Time'], axis = 1)
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])

        for x in range(0, d.shape[0]): #for each row:
            for i in range(0, d.shape[1]): #for each column:
                d.iloc[x, i] = np.round((d.iloc[x, i]/counts[x]) * total[x])

        d = d.reset_index()
        d.to_csv(self.export_path + "Takahashi2010.csv", index=False)


    def clean_sal2013(self):
        #data as downloaded from:
        #https://www.esapubs.org/archive/ecol/E094/149/

        d1 = pd.read_csv(self.import_path + "/sal2013/Table1.csv")
        d2 = pd.read_csv(self.import_path + "/sal2013/Table3.csv")

        d2.columns.values[0] = 'SampleID'

        # Create a DataFrame to hold the consolidated columns
        df_consolidated = pd.DataFrame()
        # Initialize an empty dictionary to store the new columns to be added
        new_columns_dict = {}

        # Iterate over columns and group them
        for col in d2.columns:
            base_name = col.split('.')[0]  # Extract base name (e.g., 'A' from 'A.1')
            
            # If the base_name is not already in the dictionary, calculate the sum for this group
            if base_name not in new_columns_dict:
                new_columns_dict[base_name] = d2.filter(like=base_name).sum(axis=1)

        # Convert the dictionary to a DataFrame and concatenate it with df_consolidated
        df_consolidated = pd.concat([df_consolidated, pd.DataFrame(new_columns_dict)], axis=1)

        d2 = df_consolidated
        d = pd.merge(d1, d2, on='SampleID', how='inner')
        d.columns = d.columns.str.strip()

        env_vars = ["Cruise", "Original_SampleNo", "Original_StationNo", "Daylength", "Chl", "QFChl",
        "Temperature", "QFTemp", "SurfacePAR", "QFPAR", "Kd490", "QFKd490", "PARz", "Nitrate",	
        "QFNO3", "Nitrite",	"QFNO2", "Ammonium", "QFNH4", "Phosphate", "QFPO4",	"Silicate", "QFSil", "MLD", "QFMLD"]
        d = d.drop(columns=env_vars)

        #subset only coccos from dataset:
        # Load the YAML content
        with open(self.import_path + "/sal2013/sal2013_species.yml", 'r') as file:
            data = load(file, Loader=Loader)
        # Convert to list:
        df = pd.DataFrame(list(data.items()), columns=['Species', 'Group'])
        # Create list of non-coccos:
        filtered_df = df[df['Group'] != 'coccolithophore']
        species_list = filtered_df['Species'].to_list()
        # drop any species not found in column names
        species_list = [col for col in d.columns if col in species_list]
        
        d = d.drop(columns=species_list)

        d.rename(columns = {'Lat':'Latitude'}, inplace = True)
        d.rename(columns = {'Lon':'Longitude'}, inplace = True)
        d['Method'] = "LM"

        d['Date'] = pd.to_datetime(d['Date'], dayfirst=True)
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns=['Date'], inplace=True)

        d = d[d['Latitude'].notna()]
        d = d[d['Longitude'].notna()]
        d = d[d['Depth'].notna()]
        d = d[d['Month'].notna()]
        d = d[d['Year'].notna()]
        d = d[d['Day'].notna()]
        d['Reference'] = "Sal2013"

        d = d.drop(columns=["SampleID"])

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000 #to convert from cells/ml to cells/L

        d = d.reset_index()
        d.to_csv(self.export_path + "Sal2013.csv", index=False)


    def clean_obrien2013(self):


        with open(self.yaml_path + 'references_obrien.yml', 'r') as f:
            references = load(f, Loader=Loader)

        d = pd.read_csv(self.import_path + "obrien2013.csv")

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

        d = d[d["Instrument/Method"].isin(["LM"])]

        d.rename(columns = {'Data Reference':'Reference'}, inplace = True)
        d.rename(columns = {'Depth (m)':'Depth'}, inplace = True)
        d.rename(columns = {'Longitutude':'Longitude'}, inplace = True)
        d.rename(columns = {'Instrument/Method':'Method'}, inplace = True)

        d = d.pivot_table(index=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'], columns=['Species/lowest classification'], values=['Cells/l'])
        d.columns = d.columns.droplevel(0)
        d = d.reset_index().rename_axis(None, axis=1)

        #d.reset_index(inplace=True)
        d=d.replace({"Reference": references})
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()

        d = d[d['Reference']!= 'Takahashi2000'] #wrong values
        d = d[d['Reference']!= 'meteor1929'] #method not provided for original data

        d.to_csv(self.export_path + "OBrien2013.csv", index=False)



    def clean_estrada2016(self):
        d = pd.read_csv(self.import_path + "estrada2016.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date', 'Station'], inplace=True)

        d['Method'] = "LM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv(self.export_path + "Estrada2016.csv", index=False)



    def clean_baumann2000(self):
        d = pd.read_csv(self.import_path + "Baumann2000.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).year

        d.drop(columns=['Date'], inplace=True)
        d.insert(0, 'Reference', "Baumann2000")
        d.insert(0, 'Method', 'SEM')
        d.drop_duplicates(inplace=True)
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv(self.export_path + "Baumann2000.csv", index=False)


    def clean_kleijne1984(self):
        d = pd.read_csv(self.import_path + "Kleijne1984.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d.drop(columns=['Date'], inplace=True)
        d.insert(0, 'Method', 'SEM')

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])

        d = d.reset_index()
        d.to_csv(self.export_path + "Kleijne1984.csv", index=False)


    def clean_keuter2023(self):
        d = pd.read_csv(self.import_path + "Keuter2023.csv")
        d.insert(0, 'Method', 'SEM')
        d.insert(0, 'Reference', "Keuter2023")

        d = d.set_index(['Latitude', 'Longitude', 'Depth','Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000
        d = d.reset_index()
        d.to_csv(self.export_path + "Keuter2023.csv", index=False)

    
    def clean_keuter2022(self):
        d = pd.read_csv(self.import_path + "Keuter2022.csv")
        d.insert(0, 'Method', 'SEM')
        d.insert(0, 'Reference', "Keuter2022")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d.drop(columns=['Date'], inplace=True)

        d = d.set_index(['Latitude', 'Longitude', 'Depth','Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000

        d = d.reset_index()
        d.to_csv(self.export_path + "Keuter2022.csv", index=False)

    def clean_guerreiro2023(self):
        d = pd.read_csv(self.import_path + "guerreiro2023.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date'], inplace=True)

        d['Method'] = "LM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv(self.export_path + "Guerreiro2023.csv", index=False)

    def clean_dimiza2015(self):
        d1 = pd.read_csv(self.import_path + "Dimiza2015.csv")
        d2 = pd.read_csv(self.import_path + "Dimiza2016.csv")
        d = pd.concat([d1, d2])

        d = d.assign(Reference='Dimiza2015')

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv(self.export_path + "Dimiza2015.csv", index=False)



    def preprocess_all(self):
        self.clean_estrada2016()
        self.clean_hagino200()
        self.clean_johnson2023()
        self.clean_devries2021()
        self.clean_okada1973()
        self.clean_takahashi2000()
        self.clean_cortes2001()
        self.clean_baumann2000()
        self.clean_sal2013()
        self.clean_obrien2013()
        self.clean_kleijne1984()
        self.clean_keuter2023()
        self.clean_keuter2022()
        self.clean_guerreiro2023()
        self.clean_dimiza2015()
        print("finished processing abundances")
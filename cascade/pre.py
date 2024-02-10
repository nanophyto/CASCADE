import pandas as pd
import numpy as np
from yaml import load, Loader
import glob, os
from functions import rename_synonyms, morphometric_size_estimate, volume_sphere
import math 


class pre_allometry():

    def __init__(self, import_path = "../data/unprocessed/allometry/", export_path = "../data/allometry/"):
        self.import_path = import_path
        self.export_path = export_path

    def pre_sheward_allometry(self):
        d = pd.read_csv(self.import_path + "sheward2024.csv")
        d = d.dropna()
        d = d[d['PIC pg C'] >0]
        d = d[d['Volume'] >0]
        d['species'] = d['Genus'] + " " + d['Species']
        d = d.rename(columns={"Volume": "volume", "PIC pg C": "pg pic", "Phase":"phase"})
        d = d[['species', 'volume', 'pg pic', 'phase']]
        d = rename_synonyms(d, remove_duplicate=False)
        filter = d['species'].str.contains('undefined')
        d = d[~filter]
        d.to_csv(self.export_path + "sheward2024.csv", index=False)
        print("finished processing: sheward2024")


class pre_size():

    def __init__(self, import_path = "../data/unprocessed/sizes/", export_path = "../data/sizes/"):
        self.import_path = import_path
        self.export_path = export_path


    def resample_size(self, d, species):
        d = d[d["species"]==species]

        size = []
        for i in range(len(d)):
            size_distribution = np.random.normal(d.iloc[i]['mean'], d.iloc[i]['sd'], 10000)
            size.append(size_distribution)

        species_data = {'species': [species],
                    'mean': [np.mean(size)],
                    'sd': [np.std(size)]}
        
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
        d = rename_synonyms(d)
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
        d = rename_synonyms(d)
        d['reference'] = "villiot2021b"
        d['method'] = 'literature morphometrics'
        d['n'] = 1
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d.to_csv(self.export_path + "viliot2021b.csv", index=False)

    def pre_obrien2013a(self):
        d = pd.read_csv(self.import_path + "obrien2013.csv")
        d = rename_synonyms(d)
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
        d = rename_synonyms(d)
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

        d = rename_synonyms(d, remove_duplicate=False, check_synonyms=False)

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
        d = rename_synonyms(d)
        d['reference'] = 'young2024'
        d['method'] = 'light microscopy'
        d.to_csv(self.export_path + "young2024.csv", index=False)

    def pre_gafar2019(self):
        d = pd.read_csv(self.import_path + "gafar2019.csv")
        d['mean'] = d['cell volume']
        d = rename_synonyms(d)
        d = self.resample_size_d(d)
        d['reference'] = 'gafar2019'
        d['method'] = 'light microscopy'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "gafar2019.csv", index=False)

    def pre_fiorini2011(self):
        d = pd.read_csv(self.import_path + "fiorini2011.csv")
        d = rename_synonyms(d)
        d = self.resample_size_d(d)
        d['reference'] = 'fiorini2011'
        d['method'] = 'coulter counter'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "fiorini2011.csv", index=False)

    def pre_gerecht2018(self):
        d = pd.read_csv(self.import_path + "gerecht2018.csv")
        d = rename_synonyms(d)
        d = self.resample_size_d(d)
        d['reference'] = 'gerecht2018'
        d['method'] = 'light microscopy'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "gerecht2018.csv", index=False)

    def pre_oviedo2014(self):
        d = pd.read_csv(self.import_path + "oviedo2014.csv")
        d = rename_synonyms(d)
        d = self.resample_size_d(d)
        d['reference'] = 'oviedo2014'
        d['method'] = 'coulter counter'
        d = d[['species', 'mean', 'sd', 'method', 'reference']] 
        d['n'] = 3
        d.to_csv(self.export_path + "oviedo2014.csv", index=False)

    def pre_supraha2015(self):
        d = pd.read_csv(self.import_path + "supraha2015.csv")
        d = rename_synonyms(d)
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

class pre_poc():

    def __init__(self, import_path = "../data/unprocessed/allometry/", export_path = "../data/allometry/"):
        self.import_path = import_path
        self.export_path = export_path


    def process(self):

        all_files = glob.glob(os.path.join(self.import_path, "*.csv"))

        d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
        d = d.fillna(0)

        #rename synonyms and typos:
        d = rename_synonyms(d, remove_duplicate=False, check_synonyms = False)

        return(d)




class pre_pic():

    def __init__(self, import_path = "../data/unprocessed/pic/", export_path = "../data/pic/"):
        self.import_path = import_path
        self.export_path = export_path

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
            d['pic'] = cell_pic(d['CN'], d['CL'], ks)
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
        d.to_csv(self.export_path + "sheward2014.csv", index=False)


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

        d.to_csv(self.export_path + "sheward2016.csv", index=False)


    def pre_sheward2024(self):

        #read sheward size data
        d = pd.read_csv(self.import_path + "sheward2024_pic.csv")

        d = d[d['Coccosphere condition']=="Intact"]
        d['Species'] = d['Species'].str.strip()
        d['Genus'] = d['Genus'].str.strip()
        d['species'] = d['Genus'] + " " + d['Species'] 

        d = d[['species', 'pg pic']]

        d = rename_synonyms(d, remove_duplicate=False)

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
        self.pre_devries2024()

        print("finished processing:")
        print("sheward2014")
        print("sheward2016")
        print("sheward2024")
        print("devries2024")



class merge_abundances():

    def __init__(self, import_path, export_path):
        self.export_path = export_path
        all_files = glob.glob(os.path.join(import_path, "*.csv"))

        d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)

        d = d.replace({0:pd.NA})

        d.reset_index(drop=False, inplace=True)

        d['Month'] = d['Month'].astype('int64')
        d['Year'] = d['Year'].astype('int64')

        d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'], inplace=True)

        d.rename(columns=lambda x: x.strip(), inplace=True)

        with open('/home/phyto/CoccoData/data/classification/synonyms.yml', 'r') as f:
            groupings = load(f, Loader=Loader)

        species = d.reset_index().drop(columns=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).columns

        dict = {species:k
            for k, v in groupings.items()
            for species in v}

        d = d.drop(columns=['Phaeocystis pouchetii', 'Thoracosphaera heimii'])
        d = d[d.sum(axis=1)>0]


        d = (d.rename(columns=dict)
            .groupby(level=0, axis=1, dropna=False)).sum( min_count=1).reset_index()

        d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).agg('mean')

        d = d.drop(columns=['index'])

        #check if final species are in species library
        species_library = {v: k for k, v in dict.items()}

        species_observed = d.columns

        if (set(species_observed).issubset(species_library)):
            print("all species are defined")
        else:
            raise ValueError("undefined species:" + str(set(species_observed).difference(species_library)))

        #drop any columns that contain "undefined"
        d = d[d.columns.drop(list(d.filter(regex='undefined')))]

        #drop any rows that sum to zero:
        d = d[d.sum(axis=1)>0]
        try: #make new dir if needed
            os.makedirs(self.export_path)
        except:
            None


        counts =pd.DataFrame({'count': np.count_nonzero(d.fillna(0), axis=0), 'species': d.columns})
        counts = counts[counts['species']!="Reticulofenestra sessilis"]
        filter = counts['species'].str.contains('undefined')
        counts = counts[~filter]
        counts = counts[counts['count']>0]

        non_zero_spp = counts[counts['count']>0]['species']

        d = d[non_zero_spp]

        self.counts = counts
        self.abundant = counts[counts['count']>=20]
        self.rare = counts[counts['count']<20]

        self.non_zero_spp = non_zero_spp

        self.d = d


    def export_obs_counts(self): 
        self.counts.sort_values(by=["count"], ascending=False).to_csv(self.export_path + "counts.csv", index=False )
        self.abundant.sort_values(by=["count"], ascending=False).to_csv(self.export_path + "abundant_species.csv", index=False )
        self.rare.sort_values(by=["count"], ascending=False).to_csv(self.export_path + "rare_species.csv", index=False )


    def export_ungridded_abundances(self):
        df = self.d #[self.non_zero_spp]
    
        df.to_csv(self.export_path +"ungridded_observations.csv")
        print("ungridded abundances exported to: " + self.export_path +"ungridded_observations.csv")

    def gridding(self):
        d = self.d[self.abundant['species']]
        d = d.reset_index()
        depth_bins = np.linspace(-1, 300, 62).astype(np.int64) 
        depth_labels = np.linspace(0, 300, 61).astype(np.int64) 
        d = d[d["Depth"] >= 0]
        d = d[d["Depth"] <= 301]

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
    
        d = d.reset_index()

        return(d)


    def export_gridded_abundances(self):
        d = self.gridding()

        d.to_csv(self.export_path +"gridded_observations.csv")
        
        print("gridded abundances exported to: " + self.export_path +"gridded_observations.csv")



    def print_stats(self):
        df = self.d[self.non_zero_spp]
        print("concatenated abundances: " +str(len(df)))

        t = pd.melt(df)
        t = t[t['value']>0]

        print("concatenated observations: " +str(len(t)))



        d = self.gridding()
        print("gridded abundances: " +str(len(d)))


        d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year'])
        t = pd.melt(d)
        t = t[t['value']>0]

        print("gridded observations: " +str(len(t)))

        


def pre_abundances():

    def clean_hagino200(): 
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Hagino2006.csv", na_values = "-")

        d['Latitude'] = d['Latitude'].str.rstrip('’S')
        d['Latitude'] = d['Latitude'].str.rstrip('’N')
        d['Latitude'] = d['Latitude'].str.replace('˚', '.')


        d['Longitude'] = d['Longitude'].str.rstrip('’E')
        d['Longitude'] = d['Longitude'].str.rstrip('’W')
        d['Longitude'] = d['Longitude'].str.replace('˚', '.')

        date = pd.to_datetime(d['Date'], format="%d/%m/%Y")

        d = d.groupby(by=d.columns, axis=1).sum()

        d['Reference'] = 'Hagino2006'
        d['Date'] = date
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns='Date', inplace=True)

        d = d.replace('--',np.NaN)

        d['Method'] = "SEM"
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.fillna(0, inplace=True)

        d = d.reset_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.to_csv("/home/phyto/CoccoData/data/abundances/Hagino2006.csv", index=False)



    def clean_devries2021(): 
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/deVries-etal_2020.csv")

        d['Date/Time'] = pd.to_datetime(d['Date/Time'])

        d = d[d.Reference != 'Takahashi & Okada (2000)'] #wrong data

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

        d['Method'] = "SEM"

        with open('/home/phyto/CoccoData/data/unprocessed/references/refs_devries.yml', 'r') as f:
            short_refs = load(f, Loader=Loader)

        d=d.replace({"Reference": short_refs})

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.reset_index(inplace=True)

        for i in range(len(d['Reference'].unique())):
            reference_to_export = d['Reference'].unique()[i]
            d[d.Reference==reference_to_export].to_csv("/home/phyto/CoccoData/data/abundances/" + reference_to_export + ".csv", index=False)

    def clean_okada1973(): 


        d1 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/Okada1973.csv")
        d2 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/Okada1973B.csv")
        d3 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse2.csv")
        d4 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse4.csv")
        d5 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse4_A.csv")
        d6 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse5.csv")
        d7 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse13.csv")

        d = pd.concat([d1, d2, d3, d4, d5, d6, d7])

        d['Method'] = "SEM"
        d['Reference'] = "Okada1973"
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
        d.to_csv("/home/phyto/CoccoData/data/abundances/Okada1973.csv", index=False)



    def clean_cortes2001():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/cortes2001.csv", na_values = "ND")
        
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date'], inplace=True)
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Cortes2001.csv", index=False)


    def clean_takahashi2000():

        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Takahashi_raw.csv")
        d['Reference'] = "Takahashi2000"
        d['Method'] = "SEM"
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
        d.to_csv("/home/phyto/CoccoData/data/abundances/Takahashi2010.csv", index=False)


    def clean_sal2013():

        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/sal2013.csv")
        d['Reference'] = "Sal2013"
        d.rename(columns = {'Lat':'Latitude'}, inplace = True)
        d.rename(columns = {'Lon':'Longitude'}, inplace = True)
        d['Method'] = "LM"

        d['Date'] = pd.to_datetime(d['Date'])
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

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Sal2013.csv", index=False)


    def clean_obrien2013():


        with open('/home/phyto/CoccoData/data/unprocessed/references/references_obrien.yml', 'r') as f:
            references = load(f, Loader=Loader)

        #reverse key-binding to: long-short
        #references = {v: k for k, v in references.items()}


        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/obrien2013.csv")

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

#        d = d.replace([
#                ], "unknown")
        

        print(d['Instrument/Method'].unique())

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

        d.to_csv("/home/phyto/CoccoData/data/abundances/OBrien2013.csv", index=False)



    def clean_estrada2016():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Malaspina2010.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date', 'Station'], inplace=True)

        d['Method'] = "LM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Estrada2016.csv", index=False)



    def clean_baumann2000():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Baumann2000.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).year

        d.drop(columns=['Date'], inplace=True)
        d['Reference'] = "Baumann2000"
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Baumann2000.csv", index=False)


    def clean_kleijne1984():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Kleijne1984.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d.drop(columns=['Date'], inplace=True)
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])

        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Kleijne1984.csv", index=False)


    def clean_keuter2023():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Keuter2023.csv")
        d['Method'] = "SEM"
        d['Reference'] = "Keuter2023"

        d = d.set_index(['Latitude', 'Longitude', 'Depth','Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Keuter2023.csv", index=False)

    
    def clean_keuter2022():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Keuter2022.csv")
        d['Method'] = "SEM"
        d['Reference'] = "Keuter2022"

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d.drop(columns=['Date'], inplace=True)

        d = d.set_index(['Latitude', 'Longitude', 'Depth','Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000

        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Keuter2022.csv", index=False)
    
    clean_estrada2016()
    clean_hagino200()
    clean_devries2021()
    clean_okada1973()
    clean_takahashi2000()
    clean_cortes2001()
    clean_baumann2000()
    clean_sal2013()
    clean_obrien2013()
    clean_kleijne1984()
    clean_keuter2023()
    clean_keuter2022()
import numpy as np 
import pandas as pd
import sys

sys.path.insert(0, '/home/phyto/CoccoData/coccodata/')
from functions import rename_synonyms

def pre_devries2024():
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
    d.to_csv("/home/phyto/CoccoData/data/pic/devries2024.csv", index=False)


def pre_sheward2014():

    def cell_pic(cn, cl, ks):
        pic = (2.7*ks*cn*cl**3)
        return(pic)

    df = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/pic/sheward2014.csv")

    def estimate_pic(d, species, ks):
        d = d[d['species']==species]
        d['pic'] = cell_pic(d['CN'], d['CL'], ks)
        return(d[['pic', 'species']])    

    d = pd.concat([estimate_pic(df, 'Coccolithus pelagicus', ks=0.06), 
                   estimate_pic(df, 'Coccolithus braarudii', ks=0.06)])

    d = d.groupby(by="species").agg(["mean", "std", "count"]).reset_index()

    d['sd'] = np.round(d['pic']['std'], 1)
    d['mean'] = np.round(d['pic']['mean'], 1)
    d['n'] = d['pic']['count']

    d = d[['species', 'sd', 'mean', 'n']]
    d.reset_index(inplace=True, drop=True)
    d = d.dropna()

    d['reference'] = 'sheward2014'
    d['method'] = "lab morphometrics"
    d.to_csv("/home/phyto/CoccoData/data/pic/sheward2014.csv", index=False)


def pre_sheward2016():

    def cell_pic(cn, cl, ks):
        pic = (2.7*ks*cn*cl**3)
        return(pic)

    def estimate_quadriperforatus():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/pic/sheward2017_quadriperforatus.csv")
        d['pic'] = cell_pic(d['CN'], d['CL'], 0.08)
        return(d[['pic', 'species']])

    def estimate_leptoporus():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/pic/sheward2017_leptoporus.csv")
        d['pic'] = cell_pic(d['CN'], d['CL'], 0.08)
        return(d[['pic', 'species']])

    def estimate_carteri():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/pic/sheward2017_carteri.csv")
        d['pic'] = cell_pic(d['CN'], d['CL'], 0.05)
        return(d[['pic', 'species']])

    d = pd.concat([estimate_quadriperforatus(), estimate_leptoporus(), estimate_carteri()])

    d = d.groupby(by="species").agg(["mean", "std", "count"]).reset_index()

    d['sd'] = np.round(d['pic']['std'], 1)
    d['mean'] = np.round(d['pic']['mean'], 1)
    d['n'] = d['pic']['count']
    d = d[['species', 'sd', 'mean', 'n']]

    d['reference'] = 'sheward2016'
    d['method'] = "lab morphometrics"

    d = d.dropna()
    d.reset_index(inplace=True, drop=True)
    d.to_csv("/home/phyto/CoccoData/data/pic/sheward2016.csv", index=False)


def pre_sheward2024():

    #read sheward size data
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/pic/sheward2024_pic.csv")

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
    d.to_csv("/home/phyto/CoccoData/data/pic/sheward2024.csv", index=False)


pre_sheward2024()
pre_sheward2016()
pre_sheward2014()
pre_devries2024()
print("fin")



import pandas as pd
import numpy as np
from functions import ratio_bootstrap
from math import log10, floor
from functions import rename_synonyms, bayes_bootstrap
import glob, os

class cellular_dataset_table:
    def __init__(self, d): 
        self.d = d

    def round_sig_fig(self, x, N=1):
        N = N-1
        
        return round(x, N -int(floor(log10(abs(x)))))

    def merge_results(self, species):

        df = self.d[self.d['species']==species]

        poc = df[df['variable']=='pg poc']['value']
        pic = df[df['variable']=='pg pic']['value']
        size = df[df['variable']=='diameter (um)']['value']

        size_ci_lo = np.percentile(size, 2.5)
        size_ci_up = np.percentile(size, 97.5)
        size_mean = np.percentile(size, 50)

        size_ci_lo = self.round_sig_fig(size_ci_lo, N=3)
        size_ci_up = self.round_sig_fig(size_ci_up, N=3)
        size_mean = self.round_sig_fig(size_mean, N=3)

        pic_ci_lo =  np.percentile(pic, 2.5) 
        pic_ci_up = np.percentile(pic, 97.5) 
        pic_mean = np.percentile(pic, 50) 

        pic_ci_lo = self.round_sig_fig(pic_ci_lo, N=3)
        pic_ci_up = self.round_sig_fig(pic_ci_up, N=3)
        pic_mean = self.round_sig_fig(pic_mean, N=3)


        poc_ci_lo =  np.percentile(poc, 2.5) 
        poc_ci_up = np.percentile(poc, 97.5) 
        poc_mean = np.percentile(poc, 50) 


        poc_ci_lo = self.round_sig_fig(poc_ci_lo, N=3)
        poc_ci_up = self.round_sig_fig(poc_ci_up, N=3)
        poc_mean = self.round_sig_fig(poc_mean, N=3)

        pic_poc = ratio_bootstrap(pic, poc)
        pic_poc_mean = np.percentile(pic_poc, 50) 
        pic_poc_ci_up = np.percentile(pic_poc, 97.5) 
        pic_poc_ci_lo = np.percentile(pic_poc, 2.5)
        pic_poc_ci_lo = self.round_sig_fig(pic_poc_ci_lo, N=3)
        pic_poc_ci_up = self.round_sig_fig(pic_poc_ci_up, N=3)
        pic_poc_mean = self.round_sig_fig(pic_poc_mean, N=3)   


        counts = pd.read_csv("./data/output/counts.csv")
        n_samples = counts[counts['species']==species]['count']

        d = pd.DataFrame({'species': [species], 
                        'size (med, [95% CI])':str(size_mean) + " " + str([size_ci_lo, size_ci_up]),
                        'POC (med, [95% CI])':str(poc_mean) + " " + str([poc_ci_lo, poc_ci_up]),
                        'PIC (med, [95% CI])':str(pic_mean) + " " + str([pic_ci_lo, pic_ci_up]),
                        'PIC:POC (med, [95% CI])':str(pic_poc_mean) + " " + str([pic_poc_ci_lo, pic_poc_ci_up]),


                        #, 'abundance obs':n_samples
                        })
        

        return(d)

    def concat_estimates(self):

        df = []
        species_list = self.d['species'].unique()
        for i in range(len(species_list)):
            print(species_list[i])
            df.append(self.merge_results(species_list[i]))
            print("finished appending species #" + str(i+1) + 
                  " out of " + str(len(species_list)) + " species")

        df = pd.concat(df)

        df = df.sort_values(by=['species'])
        
        print(df.to_latex(index=False,
                        formatters={"name": str.upper},
                        float_format="{:.1f}".format,
        )) 


class gridded_datasets():

    def __init__(self, abundances_path, carbon_path, export_path):

        self.d = pd.read_csv(abundances_path)
        self.df = pd.read_csv(carbon_path)

        #then merge with abundances
        species_list = self.d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year']).columns

        gridded_poc = []
        gridded_pic = []

        for i in range(0, len(species_list)):
        #    print(species_list[i])
            gridded_poc.append(self.merge_abundance_c(species_list[i], variable = 'pg poc'))
            gridded_pic.append(self.merge_abundance_c(species_list[i], variable = 'pg pic'))
        
        gridded_poc = pd.concat(gridded_poc)
        gridded_pic = pd.concat(gridded_pic)

        gridded_poc.to_csv(export_path + "gridded_poc.csv", index=False)
        gridded_pic.to_csv(export_path + "gridded_pic.csv", index=False)



    def merge_abundance_c(self, species, variable = 'pg poc'):
        df = self.df[self.df['species']==species]
        d = self.d.dropna(subset=[species])
        poc = df[df['variable']==variable]['value']

        d_poc = pd.DataFrame({
            "ci_95_lo" + " (" + variable + " L-1)" :  d[species].values*np.percentile(poc, 2.5), 
            "ci_95_up" + " (" + variable + " L-1)"  : d[species].values*np.percentile(poc, 97.5), 
            "ci_99_lo" + " (" + variable + " L-1)" :  d[species].values*np.percentile(poc, 0.5), 
            "ci_99_up" + " (" + variable + " L-1)"  : d[species].values*np.percentile(poc, 99.5), 
            "med" : d[species].values*np.percentile(poc, 50),
            "species":[species]*len(d),
            "lat": d['Latitude'],
            "lon": d['Longitude'],
            "depth": d['Depth'],
            "month": d['Month'],
            "year": d['Year']
            })

        return(d_poc)



def split_resampled_cellular_dataset(data_path = "./data/output/cellular_dataset.csv", 
                                     export_path = "./data/Zenodo/resampled_cellular_dataset/"):
    d = pd.read_csv(data_path)
    print(d)
    diameter = d[d['variable']=='diameter (um)']
    print(diameter)
    volume = d[d['variable']=='volume (um3)']
    pic = d[d['variable']=='pg pic']
    poc = d[d['variable']=='pg poc']

    pic = pic.replace('pg pic', 'cellular carbon (pg pic)', regex=True)
    poc = poc.replace('pg poc', 'cellular carbon (pg poc)', regex=True)

    pic.to_csv(export_path + "resampled_PIC.csv", index=False)
    poc.to_csv(export_path + "resampled_POC.csv", index=False)
    diameter.to_csv(export_path + "resampled_diameter.csv", index=False)
    volume.to_csv(export_path + "resampled_volume.csv", index=False)


def merge_literature_datasets(pic_path, poc_path, size_path, export_path):

    def import_data(path):

        all_files = glob.glob(os.path.join(path, "*.csv"))

        d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
        d = d.fillna(0)

        #rename synonyms and typos:
        d = rename_synonyms(d, remove_duplicate=False, check_synonyms = True)

        return(d)
    
    pic = import_data(pic_path)
    poc = import_data(poc_path)
    size = import_data(size_path)

    pic.to_csv(export_path + "literature_PIC.csv", index=False)
    poc.to_csv(export_path + "literature_POC.csv", index=False)
    size.to_csv(export_path + "literature_size.csv", index=False)

    print("finished merging literature datasets")
    print("exported to: " + export_path)

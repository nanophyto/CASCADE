import pandas as pd
import numpy as np
from functions import ratio_bootstrap
from math import log10, floor


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
        size = df[df['variable']=='diameter']['value']

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


        counts = pd.read_csv("./data/output/species_list_full.csv")
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


class merge_all():

    def __init__(self, abundances_path, carbon_path):

        self.d = pd.read_csv(abundances_path)
        self.df = pd.read_csv(carbon_path)

        #then merge with abundances
        species_list = self.d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Month', 'Year']).columns

        estimates = []

        for i in range(0, len(species_list)):
            print(species_list[i])
            estimates.append(self.merge_abundance_c(species_list[i], variable = 'pg poc'))
            estimates.append(self.merge_abundance_c(species_list[i], variable = 'pg pic'))
        
        self.estimates = pd.concat(estimates)

    def merge_abundance_c(self, species, variable = 'pg poc'):
        df = self.df[self.df['species']==species]
        d = self.d.dropna(subset=[species])
        poc = df[df['variable']==variable]['value']

        d_poc = pd.DataFrame({
            "ci_95_lo" :  d[species].values*np.percentile(poc, 2.5), 
            "ci_95_up" : d[species].values*np.percentile(poc, 97.5), 
            "med" : d[species].values*np.percentile(poc, 50),
            "species":[species]*len(d),
            "lat": d['Latitude'],
            "lon": d['Longitude'],
            "depth": d['Depth'],
            "month": d['Month'],
            "year": d['Year']
            })

        t = pd.melt(d_poc, id_vars=['species', 'lat', 'lon', 'depth', 'month', 'year'], 
                    value_vars=['ci_95_lo', 'ci_95_up', 'med'])

        t = t.sort_values(by=['month', 'year',  'depth', 'lat', 'lon'])
        return(d_poc)

    def export_csv(self, export_path):
        self.estimates.to_csv(export_path, index=False)
        print("finished exporting to:")


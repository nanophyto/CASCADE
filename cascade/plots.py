import pandas as pd
import numpy as np
import sys
import geopandas as gpd
from shapely.geometry import Point # Point class
from yaml import load, Loader
import glob, os
import pandas as pd
import numpy as np
from functions import ratio_bootstrap
from math import log10, floor
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.graphics.api import abline_plot
import cartopy.crs as ccrs
import cartopy as cart
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from main import library
from functions import rename_synonyms, bayes_bootstrap
from  matplotlib.ticker import FuncFormatter
import matplotlib.gridspec as gridspec


class LM_SEM_size:
    def __init__(self, d): 
        self.d = d

    def make_plot(self):
        d = pd.read_csv("/home/phyto/CoccoData/data/SEM_LM_counts/Bollman2002.csv")

        # Melt the DataFrame
        d = pd.melt(d, id_vars=['method', 'sample'], var_name='species', value_name='abundance')

        # Rename synonyms and take the sum
        d = rename_synonyms(d, index=['species', 'method', 'sample'], take_sum=True)

        # Filter and rename columns for 'LM'
        d1 = d[d['method']=='LM'].drop(columns=['method']).rename(columns={"abundance": "LM"})

        # Filter and rename columns for 'SEM'
        d2 = d[d['method']=='SEM'].drop(columns=['method']).rename(columns={"abundance": "SEM"})

        # Merge DataFrames on 'species' and 'sample'
        d = pd.merge(d1, d2, on=["species", "sample"]).dropna()

        d.set_index(['sample', 'species'], inplace = True)

        d = d[(d['SEM']>0) & (d['LM']>0)]

        d['ratio'] = (d['LM']-d['SEM'])/(d['SEM'])*100

        x = bayes_bootstrap(d['ratio'])
        sns.histplot(x)

        plt.title("LM vs SEM abundances (mean % difference)")
        plt.xlabel("mean % difference")

        plt.show()

        print(np.round(np.percentile(x, 2.5) , 2))
        print(np.round(np.percentile(x, 97.5) , 2))

        print(np.round(np.percentile(x, 50) , 2))

        print("fin")



class LM_SEM_size:
    def __init__(self, d): 
        self.d = d


    def make_plot(self):

        n = 10000
        d = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")
        species_list = d['species']

        m = library('/home/phyto/CoccoData/data/classification/phases.yml',
                    '/home/phyto/CoccoData/data/classification/family.yml',
                    "/home/phyto/CoccoData/data/sizes/",
                    "/home/phyto/CoccoData/data/pic/",
                    "/home/phyto/CoccoData/data/poc/",
                    species_list)

        #create a list of HOL species for future use:
        HOL_list = m.return_HOL_list()

        #return the named tuple so we can use it:
        ntpl = m.return_ntpl()

        lm_library_mean =  pd.Series([t.mean for t in ntpl  if t.id == "young2024" or t.id == "viliot2021a"])
        lm_library_spp =  pd.Series([t.species for t in ntpl  if t.id == "young2024" or t.id == "viliot2021a"])
        lm = pd.DataFrame({'species': lm_library_spp, 'mean': lm_library_mean})
        lm = lm.groupby(['species']).mean()
        lm['method'] = 'lm'

        sem_library_mean =  [t.mean for t in ntpl  if t.id == "sheward2024" or t.id == "viliot2021b"]
        sem_library_spp =  [t.species for t in ntpl  if t.id == "sheward2024"or t.id == "viliot2021b"]
        sem = pd.DataFrame({'species': sem_library_spp, 'mean': sem_library_mean})
        sem = sem.groupby(['species']).mean()
        sem['method'] = 'sem'

        d = pd.concat([sem, lm])
        d.reset_index(inplace=True)

        d = d.pivot(index='species', columns='method', values='mean')

        d['rLM'] = d['lm']/(d['lm']+d['sem'])
        d['rSEM'] = d['sem']/(d['lm']+d['sem'])



        fig, ax = plt.subplots()

        d = d.dropna()
        no_drop = len(d)
        x = np.log(d['lm'])
        y = np.log(d['sem'])
        sns.scatterplot(x=x, y=y, ax=ax, color="firebrick", s=100)
        abline_plot(slope=1, intercept=0, ax=ax, color='black', linestyle='dashed')

        ax.set_xlim(2, np.nanmax([x*1.1, y*1.1]))
        ax.set_ylim(2, np.nanmax([x*1.1, y*1.1]))

        plt.show()



class multipage_plots():

    def __init__(self, ds, species_list, n_page, out_path):
        self.species_list = species_list
        self.n = n_page
        self.out_path = out_path
        self.ds = ds

    def latlon_plots(self, species, ax=None, title=None, log_scale=False, add_colorbar=False):

        d = self.ds[self.ds[species]>0]

        d = d.groupby(['Latitude', 'Longitude'])[species].agg('mean').reset_index()

        ax.coastlines()
        ax.add_feature(cart.feature.LAND, zorder=-1,facecolor="gray")
        ax.gridlines(draw_labels=True, zorder=-1)
        ax.set_facecolor('grey')

        grid = sns.scatterplot(x=d['Longitude'], y=d['Latitude'], ax=ax, hue=np.log(d[species]), 
                               palette='viridis', edgecolor = "none", s=100, alpha=1, zorder=10) # marker="$\circ$",  ec="face", alpha=0.5
       
        ax.get_legend().remove()    
        ax.set_ylim(-90, 90)
        ax.set_xlim(-180, 180)

        if title==None:
            ax.set_title(species)
        else:
            ax.set_title(title)

    def depthtime_plots(self, species, ax=None, title=None, log_scale=False, add_colorbar=False):        
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        d = self.ds[self.ds[species]>0]
        d = d.groupby(['Depth', 'months since winter solstice'])[species].agg('mean').reset_index()

        grid = sns.scatterplot(x=d['months since winter solstice'].astype(int), 
                               y=d['Depth'], palette='viridis', 
                               hue=np.log10(d[species]), ax=ax, 
                               edgecolor = "none") 
        sm = plt.cm.ScalarMappable(cmap="viridis")

        ax.get_legend().remove()
        ax.figure.colorbar(sm)
        ax.set_facecolor('grey')
        ax.set_ylim(300, 0)
        ax.set_xlim(1, 12)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))

        if title==None:
            ax.set_title(species)
        else:
            ax.set_title(title)


    def sub_plots_ax(self, species):
        gs = gridspec.GridSpec(6,3, width_ratios=[1,2,2], 
                               height_ratios=[1, 1, 1, 1, 1, 1])

        cols = 3
        rows = len(species)
        plt.figure(figsize=(18, 28))
        species = list(species)
        for j in range(rows):

            ax1 = subplot(rows,cols,(j*cols)+1)
            ax2 = subplot(rows,cols,(j*cols)+2, projection=ccrs.PlateCarree())
            ax3 = subplot(rows,cols,(j*cols)+3)

            self.latlon_plots(ax=ax2, species=species[j])
            self.depthtime_plots(ax=ax3, species=species[j])

            try:
                image = plt.imread("/home/phyto/CASCADE/data/sem/" + species[j] + ".jpeg")
            except:
                None
            ax1.imshow(image)  
            ax1.set_axis_off()
            ax2.axes.get_yaxis().set_visible(False)

            species_title =  str(species[j]) 
            ax1.title.set_text(species_title)

            ax2_title = 'Latitude vs Longitude'
            ax2.title.set_text(ax2_title )
            
            ax3_title = 'Time vs Depth'
            ax3.title.set_text(ax3_title)
            
        plt.subplots_adjust(left=None, bottom=None, right=None, 
                            top=None, wspace=0.5, hspace=1)
        
    def export_pdf(self):

        n = self.n
        final = [self.species_list[i * n:(i + 1) * n] for i in range((len(self.species_list) + n - 1) // n )]  

        with PdfPages(self.out_path) as pdf:
            for i in range(len(final)):
                print(i)
                self.sub_plots_ax(final[i])
                pdf.savefig(pad_inches=0.02, bbox_inches="tight")
                plt.close()

        print("finished exporting pdf to: " + self.out_path)




class sd_simulations:

    def __init__(self):
        None

    def sim_sd_minmax(self, k):

        def sd_minmax(n=100, generator=np.random.normal):
            sample = generator(size=n)
            real = np.std(sample)
            minmax = (sample.max() - sample.min())/k
            return (real, minmax, ((minmax - real)*100)/real)
            
        all_sims = []
        for n in np.arange(1, 30, 1):
            sims = [sd_minmax(n=n) for _ in range(1000)]
            all_sims.append(np.row_stack(sims))

        all_simdf = pd.concat([
            pd.DataFrame(sim).assign(n=n) 
            for (sim,n) in 
            zip(all_sims, np.arange(1, 30, 1))
            ])
        all_simdf.columns = ['real','minmax','% error', 'n_samples']
        
        all_simdf['k'] = k 
        return(all_simdf)
    
    def plot(self):

        all_simdf = pd.concat([self.sim_sd_minmax(k=2), self.sim_sd_minmax(k=4)])

        g = sns.FacetGrid(all_simdf, col="k")

        g.map_dataframe(sns.boxplot, x='n_samples', y='% error')

        ax1, ax2 = g.axes[0]

        ax1.axhline(0, ls='--')
        ax2.axhline(0, ls='--')

        ax1.text(0.05, 0.95, "A)", transform=ax1.transAxes,
            fontsize=16, fontweight='bold', va='top')
        ax2.text(0.05, 0.95, "B)", transform=ax2.transAxes,
            fontsize=16, fontweight='bold', va='top')

        ax1.tick_params(axis='x', labelrotation=90)
        ax2.tick_params(axis='x', labelrotation=90)

        #f.autofmt_xdate(rotation=90)
        #ax.legend(ncol=2, title='Sample Size')
        plt.show()


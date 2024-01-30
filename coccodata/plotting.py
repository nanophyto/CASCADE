import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import scipy.stats as st 
import cartopy.crs as ccrs
import cartopy as cart
from matplotlib import cm


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

class multipage_plots():

    def __init__(self, ds, species_list, n_page, out_path):
        self.species_list = species_list
        self.n = n_page
        self.out_path = out_path
        self.ds = ds

    def latlon_plots(self, species, ax=None, title=None, log_scale=False, add_colorbar=False):

        #ax = plt.axes(projection=ccrs.PlateCarree())
        d = self.ds[self.ds[species]>0]

        d = d.groupby(['Latitude', 'Longitude'])[species].agg('mean').reset_index()

        ax.coastlines()
        ax.add_feature(cart.feature.LAND, zorder=1, edgecolor='k', facecolor="gray")
        ax.gridlines(draw_labels=True, zorder=1)
        

        grid = sns.scatterplot(x=d['Longitude'], y=d['Latitude'], ax=ax, hue=np.log(d[species]), 
                               palette='viridis', edgecolor = "none") # marker="$\circ$",  ec="face", alpha=0.5
       
        #sm = plt.cm.ScalarMappable(cmap="viridis")
        ax.get_legend().remove()    
        #ax.figure.colorbar(sm)

        # plt.scatter(, 
        #             s=100,
        #             alpha=0.75,
        #             transform=ccrs.PlateCarree()) ## Important

        ax.set_ylim(-90, 90)
        ax.set_xlim(-180, 180)


        #cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])#
        if title==None:
            ax.set_title(species)
        else:
            ax.set_title(title)

    def depthtime_plots(self, species, ax=None, title=None, log_scale=False, add_colorbar=False):

        #ax = plt.axes(projection=ccrs.PlateCarree())
        
        d = self.ds[self.ds[species]>0]

        d = d.groupby(['Depth', 'months since winter solstice'])[species].agg('mean').reset_index()


        grid = sns.scatterplot(x=d['months since winter solstice'], y=d['Depth'], palette='viridis', 
                               hue=np.log10(d[species]), ax=ax, edgecolor = "none") #marker="$\circ$",  ec="face", 
        sm = plt.cm.ScalarMappable(cmap="viridis")
        ax.get_legend().remove()
        
    
        ax.figure.colorbar(sm)


        ax.set_ylim(300, 0)
        ax.set_xlim(1, 12)

        #cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])#
        if title==None:
            ax.set_title(species)
        else:
            ax.set_title(title)


    def sub_plots_ax(self, species):
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
                image = plt.imread("/home/phyto/CoccoData/data/sem/" + species[j] + ".jpeg")
            except:
                image = plt.imread("/home/phyto/SEM_eva/JRYSEM-287-07p.jpeg")
            ax1.imshow(image)  
            ax1.set_axis_off()
            ax2.axes.get_yaxis().set_visible(False)


            species_title =  str(species[j]) 
            ax1.title.set_text(species_title)

            ax2_title = r'\textbf{Diameter (um)'
            ax2.title.set_text(ax2_title )
            
            ax3_title = r'\textbf{Particulate Carbon (pg C)}'
            ax3.title.set_text(ax3_title)
            

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)
    


    def export_pdf(self):

        n = self.n
        final = [self.species_list[i * n:(i + 1) * n] for i in range((len(self.species_list) + n - 1) // n )]  

        with PdfPages(self.out_path) as pdf:
            for i in range(len(final)):
                print(i)
                self.sub_plots_ax(final[i])
                pdf.savefig()
                plt.close()

        print("finished exporting pdf to: " + self.out_path)


d = pd.read_csv("/home/phyto/CoccoData/gridded_abundances.csv")
d['months since winter solstice'] = d['Month']
d.loc[(d['Latitude'] <0) & (d['Month'] == 1),'months since winter solstice']= 7
d.loc[(d['Latitude'] <0) & (d['Month'] == 2),'months since winter solstice']= 8
d.loc[(d['Latitude'] <0) & (d['Month'] == 3),'months since winter solstice']= 9
d.loc[(d['Latitude'] <0) & (d['Month'] == 4),'months since winter solstice']= 10
d.loc[(d['Latitude'] <0) & (d['Month'] == 5),'months since winter solstice']= 11
d.loc[(d['Latitude'] <0) & (d['Month'] == 6),'months since winter solstice']= 12
d.loc[(d['Latitude'] <0) & (d['Month'] == 7),'months since winter solstice']= 1
d.loc[(d['Latitude'] <0) & (d['Month'] == 8),'months since winter solstice']= 2
d.loc[(d['Latitude'] <0) & (d['Month'] == 9),'months since winter solstice']= 3
d.loc[(d['Latitude'] <0) & (d['Month'] == 10),'months since winter solstice']= 4
d.loc[(d['Latitude'] <0) & (d['Month'] == 11),'months since winter solstice']= 5
d.loc[(d['Latitude'] <0) & (d['Month'] == 12),'months since winter solstice']= 6


ds = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")
ds = ds.sort_values(by=['species'])
#species_list = list(ds['species'])[0:6]
species_list = list(ds['species'])

# d = pd.read_csv("/home/phyto/CoccoData/final_pic.csv")

# # #m = multipage_plots(d=d, species_list=species_list[0:7], n_page=8, out_path='/home/phyto/library.pdf')

m = multipage_plots(ds=d, species_list=species_list, n_page=6, out_path='/home/phyto/library.pdf')
m.export_pdf()

# print("fin")

# species_list = ["Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus"]

# multipage_plots(d=d, species_list=species_list, n_page=8, out_path='/home/phyto/multipage_pdf.pdf')

#print("fin")
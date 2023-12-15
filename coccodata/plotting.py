import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from pylab import *
import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

class multipage_plots():

    def __init__(self, d, species_list, n_page, out_path):
        self.d = d
        self.species_list = species_list
        self.n = n_page
        self.out_path = out_path

    def dist_plots(self, ax, variable, species, log_scale=False):
        df = self.d[self.d['species']==species]
        poc = df[df['variable']=='pg poc']['value']
        pic = df[df['variable']=='pg pic']['value']
        df = df[df['variable']==variable]

        
        try:
            # if percentile is not None:
            #     d = d[d['value']> np.percentile(df['value'], 100-percentile)]
            #     d = d[d['value']< np.percentile(df['value'], percentile)]
            x = np.random.choice(df['value'], 1000)
            #x = d['value']
        except:
            x = []

        #grid = sns.histplot(x = x,ax=ax, kde = True, stat="density")
        grid = sns.kdeplot(x = x, ax=ax, log_scale=log_scale)


        if variable == 'pg poc':
            #pic_poc = np.random.choice(pic, 1000)/np.random.choice(poc, 1000)


            pic_median = np.round(np.median(pic), 1)
            poc_median = np.round(np.median(poc), 1)
            #pic_poc_median = np.round(np.median(pic_poc), 1)

            pic_ci = np.round(np.std(pic), 1)
            poc_ci = np.round(np.std(poc), 1)
            #pic_poc_ci = np.round(np.std(pic_poc), 1)

            # Build a rectangle in axes coords

            grid.text(0.95, 0.95, r'POC:' +  str(poc_median) + '± ' + str(poc_ci),
                horizontalalignment='right',
                verticalalignment='top',
                transform=grid.transAxes)
            grid.text(0.95, 0.85, r'PIC:' +  str(pic_median) + '± ' + str(pic_ci),
                horizontalalignment='right',
                verticalalignment='top',
                transform=grid.transAxes)        
            #grid.text(0.95, 0.75, r'PIC/POC:' +  str(pic_poc_median) + '± ' + str(pic_poc_ci),
            #    horizontalalignment='right',
            #    verticalalignment='top',
            #    transform=grid.transAxes)   
            
            p = plt.Rectangle((0, 0), 1, 1, fill=False)
            p.set_transform(grid.transAxes)
            p.set_clip_on(False)
            grid.add_patch(p)
        else:
            None
        return(grid)

    def sub_plots_ax(self, species):
        cols = 4
        rows = len(species)
        plt.figure(figsize=(18, 28))
        species = list(species)
        for j in range(rows):

            ax1 = subplot(rows,cols,(j*cols)+1)
            ax2 = subplot(rows,cols,(j*cols)+2)
            ax3 = subplot(rows,cols,(j*cols)+3)
            ax4 = subplot(rows,cols,(j*cols)+4)

            self.dist_plots(ax=ax2, variable="diameter", species=species[j])
            self.dist_plots(ax=ax3, variable="pg poc", species=species[j], log_scale=True)
            self.dist_plots(ax=ax3, variable="pg pic", species=species[j], log_scale=True)
            self.dist_plots(ax=ax4, variable="pg pic", species=species[j], log_scale=True)
            try:
                image = plt.imread("/home/phyto/CoccoData/data/sem/" + species[j] + ".jpeg")
            except:
                image = plt.imread("/home/phyto/SEM_eva/JRYSEM-287-07p.jpeg")
            ax1.imshow(image)  
            ax1.set_axis_off()
            ax2.axes.get_yaxis().set_visible(False)
            ax3.axes.get_yaxis().set_visible(False)
            ax4.axes.get_yaxis().set_visible(False)
            species_title = r'\textbf{\textit{' + str(species[j]) + '}}'
            ax1.title.set_text(species_title)
            ax2_title = r'\textbf{Diameter (um)'

            ax2.title.set_text(ax2_title )
            
            ax3_title = r'\textbf{Particulate Carbon (pg C)}'
            ax3.title.set_text(ax3_title)
            
            ax4.title.set_text("PIC (pg C)")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)
        


    def export_pdf(self):

        n = self.n
        final = [species_list[i * n:(i + 1) * n] for i in range((len(self.species_list) + n - 1) // n )]  

        with PdfPages(self.out_path) as pdf:
            for i in range(len(final)):
                print(i)
                self.sub_plots_ax(final[i])
                pdf.savefig()
                plt.close()

        print("finished exporting pdf to: " + self.out_path)




d = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")
species_list = d['species']
d= pd.read_csv("/home/phyto/CoccoData/fml2.csv")
#m = multipage_plots(d=d, species_list=species_list[0:7], n_page=8, out_path='/home/phyto/library.pdf')
m = multipage_plots(d=d, species_list=species_list, n_page=6, out_path='/home/phyto/library.pdf')

m.export_pdf()

print("fin")#

# species_list = ["Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus"]

# multipage_plots(d=d, species_list=species_list, n_page=8, out_path='/home/phyto/multipage_pdf.pdf')

#print("fin")
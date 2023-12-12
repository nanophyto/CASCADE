import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from pylab import *
import numpy as np


def multipage_plots(d, species_list, n_page, out_path):

    def dist_plots(ax, variable, species):
        df = d[d['Species']==species]
        grid = sns.histplot(df, x=variable, ax=ax)
        return(grid)

    def sub_plots_ax(species):
        cols = 4
        rows = len(species)

        for j in range(rows):

            ax1 = subplot(rows,cols,(j*cols)+1)
            ax2 = subplot(rows,cols,(j*cols)+2)
            ax3 = subplot(rows,cols,(j*cols)+3)
            ax4 = subplot(rows,cols,(j*cols)+4)

            dist_plots(ax=ax1, variable="diameter", species=species[j])
            dist_plots(ax=ax2, variable="volume", species=species[j])
            dist_plots(ax=ax3, variable="pg poc", species=species[j])
            image = plt.imread("/home/phyto/SEM_eva/JRYSEM-287-07p.jpeg")
            ax4.imshow(image)    



    n = n_page
    # using list comprehension 
    final = [species_list[i * n:(i + 1) * n] for i in range((len(species_list) + n - 1) // n )]  

    with PdfPages(out_path) as pdf:
        for i in range(len(final)):
            sub_plots_ax(final[i])
            pdf.savefig()

    print("finished exporting pdf to: " + out_path)


d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
species_list = ["Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra"]

multipage_plots(d=d, species_list=species_list, n_page=6, out_path='/home/phyto/multipage_pdf.pdf')

print("fin")
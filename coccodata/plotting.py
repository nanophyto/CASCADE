import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from pylab import *
import numpy as np


def multipage_plots(d, species_list, n_page, out_path):

    def dist_plots(ax, variable, species):
        df = d[d['species']==species]
        df = df[df['variable']==variable]
        grid = sns.histplot(df, x="value", ax=ax)
        return(grid)

    def sub_plots_ax(species):
        cols = 4
        rows = len(species)
        plt.figure(figsize=(18, 28))

        for j in range(rows):

            ax1 = subplot(rows,cols,(j*cols)+1)
            ax2 = subplot(rows,cols,(j*cols)+2)
            ax3 = subplot(rows,cols,(j*cols)+3)
            ax4 = subplot(rows,cols,(j*cols)+4)

            dist_plots(ax=ax2, variable="diameter", species=species[j])
            dist_plots(ax=ax3, variable="volume", species=species[j])
            dist_plots(ax=ax4, variable="pg poc", species=species[j])
            image = plt.imread("/home/phyto/SEM_eva/JRYSEM-287-07p.jpeg")
            ax1.imshow(image)  
            ax1.set_axis_off()
            ax2.axes.get_yaxis().set_visible(False)
            ax3.axes.get_yaxis().set_visible(False)
            ax4.axes.get_yaxis().set_visible(False)

            ax1.title.set_text(species[j])
            ax2.title.set_text("diameter (um)")
            ax3.title.set_text("volume (um^3)")
            ax4.title.set_text("POC (pg C)")

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
        




    n = n_page
    final = [species_list[i * n:(i + 1) * n] for i in range((len(species_list) + n - 1) // n )]  

    with PdfPages(out_path) as pdf:
        for i in range(len(final)):
            sub_plots_ax(final[i])
            pdf.savefig()
            plt.close()

    print("finished exporting pdf to: " + out_path)

#
# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
# species_list = ["Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus",
#                 "Emiliania huxleyi", "Calcidiscus leptoporus", "Syracosphaera pulchra", "Calcidiscus leptoporus"]

# multipage_plots(d=d, species_list=species_list, n_page=8, out_path='/home/phyto/multipage_pdf.pdf')

#print("fin")
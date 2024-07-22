import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.graphics.api import abline_plot
import cartopy.crs as ccrs
import cartopy as cart
from pylab import *
from main import library
from functions import rename_synonyms, bayes_bootstrap, months_since_winter_solstice_df
import matplotlib as mpl


def LM_SEM_size_plot(d_path):

    d = pd.read_csv(d_path)

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

    fig, ax = plt.subplots()

    x = bayes_bootstrap(d['ratio'])
    sns.histplot(x, ax=ax)

    plt.title("LM vs SEM abundances (mean % difference)")
    plt.xlabel("mean % difference")

    plt.show()

    print(np.round(np.percentile(x, 2.5) , 2))
    print(np.round(np.percentile(x, 97.5) , 2))

    print(np.round(np.percentile(x, 50) , 2))

    print("fin")
    return(fig)



def LM_SEM_size(log=True, figsize=(12, 12)):

    d = pd.read_csv("./data/output/counts.csv")
    species_list = d['species']

    m = library('./data/classification/phases.yml',
                './data/classification/family.yml',
                "./data/sizes/",
                "./data/pic/",
                "./data/poc/",
                species_list)

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


    d['ratio'] = d['sem']/d['lm']
    d.reset_index(inplace=True)
    d_copy = d.copy()

    d1 = d_copy[d_copy['ratio']<=0.2]
    d2 = d_copy[d_copy['ratio']>=5]
    d_dropped = pd.concat([d1, d2])

    d = d[d['ratio']>0.2]
    d = d[d['ratio']<5]


    if log==False:
        x = d['lm']
        y = d['sem']
        x_dropped = d_dropped['lm']
        y_dropped = d_dropped['sem']
    else:
        x = np.log10(d['lm'])
        y = np.log10(d['sem'])
        x_dropped = np.log10(d_dropped['lm'])
        y_dropped = np.log10(d_dropped['sem'])


    fig, ax = plt.subplots(figsize=figsize)

    abline_plot(slope=1, intercept=0, ax=ax, color='black', linestyle='dashed')
    sns.scatterplot(x=x, y=y, ax=ax, s=100)
    sns.scatterplot(x=x_dropped, y=y_dropped, ax=ax, s=100, color="red")

    if log==False:
        ax.set_xlim(0, np.nanmax([x, y])+100)
        ax.set_ylim(0, np.nanmax([x, y])+100)
    else:
        ax.set_xlim(0.5, np.nanmax([x, y])+1)
        ax.set_ylim(0.5, np.nanmax([x, y])+1)

    fig.suptitle('Cell size estimates SEM vs LM', weight='bold')

    plt.ylabel("SEM size (um, log10)")
    plt.xlabel("LM size (um, log10)")

    return(fig)
#    plt.show()



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
            ], ignore_index=True)
        all_simdf.columns = ['real','minmax','% error', 'n_samples']
        
        all_simdf['k'] = k 
        return(all_simdf)
    
    def plot(self):

        all_simdf = pd.concat([self.sim_sd_minmax(k=2), self.sim_sd_minmax(k=4)])

#        g = sns.FacetGrid(all_simdf, col="k")
        g = sns.FacetGrid(all_simdf, col="k", height=5, aspect=1.2)

        g.map_dataframe(sns.boxplot, x='n_samples', y='% error')

        ax1, ax2 = g.axes[0]

        ax1.axhline(0, ls='--')
        ax2.axhline(0, ls='--')

        ax1.text(0.05, 0.95, "a)", transform=ax1.transAxes,
            fontsize=16, fontweight='bold', va='top')
        ax2.text(0.05, 0.95, "b)", transform=ax2.transAxes,
            fontsize=16, fontweight='bold', va='top')

        ax1.tick_params(axis='x', labelrotation=90)
        ax2.tick_params(axis='x', labelrotation=90)

        #f.autofmt_xdate(rotation=90)
        #ax.legend(ncol=2, title='Sample Size')
        return(g)


def depth_time_samples_plot(d):
    df = months_since_winter_solstice_df(d)

    # Bin the data
    depth_bins = np.arange(0, 305, 25)  # Bins for depth (0, 25, 50, 75, 100)
    month_bins = np.arange(0, 13, 1)  # Bins for months (0, 1, 2, ..., 12)

    # Bin the data
    #df['Depth'] = pd.cut(df['Depth'], bins=depth_bins, right=False)
    #df['months_since_winter_solstice'] = pd.cut(df['months_since_winter_solstice'], bins=month_bins, right=False)
    df = df.sort_values(
        by="months_since_winter_solstice")

    df['Depth'] = pd.cut(df['Depth'], bins=depth_bins, right=False, labels=[f'{depth_bins[i]}-{depth_bins[i+1]}' for i in range(len(depth_bins)-1)])
    #df['months_since_winter_solstice'] = pd.cut(df['months_since_winter_solstice'], bins=month_bins, right=False, labels=[f'{month_bins[i]}-{month_bins[i+1]}' for i in range(len(month_bins)-1)])
    df['months_since_winter_solstice']  = df['months_since_winter_solstice'].astype(str)

    cmap = mpl.cm.Blues(np.linspace(0,1,20))
    darkened_blues = mpl.colors.ListedColormap(cmap[5:,:-1])
    # Create the joint plot
    g = sns.jointplot(
        x='months_since_winter_solstice',
        y='Depth',
        data=df,
        ratio=10,
        #height=25,
        kind='hist',  # Use hist to get both 2D histogram and marginal histograms
        marginal_kws={'bins': month_bins, 'color': 'steelblue'},  # Bins for marginal histograms
        cmap=darkened_blues  # Color map for the 2D histogram
    )

    # Add labels and title
    g.set_axis_labels('Months Since Winter Solstice', 'Depth (m)')

    g.fig.set_figwidth(24)
    g.fig.set_figheight(8)
    #g.savefig("filename.png", dpi=300)
    return(g)
    # Show the plot
#    plt.show()
#    plt.rcParams.update({'font.size':20})



def plot_samples_latlon(d, fill=None, ax=None, fig=None, log=False, title=""):

    if fig==None:
        fig = plt.figure(figsize=(20, 10))
    if ax==None:
        projection = ccrs.Robinson(central_longitude=-160)
        ax= plt.axes(projection=projection)

    ax.coastlines()
    ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k', facecolor="lightgrey")
    ax.gridlines(draw_labels=True,)

    if fill==None:
        color="dodgerblue"
    else:
        if log==False:
            color=d[fill]
        else:
            color=np.log(d[fill])


    plt.scatter(x=d['Longitude'], y=d['Latitude'],
                c=color,
                s=100,
                alpha=0.75,
                transform=ccrs.PlateCarree()) ## Important

#    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])#
    ax.set_title(title)
#    plt.colorbar(p, cax=cax)

import pandas as pandas
import numpy as np 
from collections import namedtuple
from yaml import load, Loader
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.insert(0, '/home/phyto/CoccoData/coccodata/')
from library import library
#from regression import regression_simulation
#from plotting import multipage_plots
from group import estimate_group_pic_poc
from statsmodels.graphics.api import abline_plot
from scipy import stats
plt.rcParams.update({'font.size': 18})

n = 10000

d = pd.read_csv("/home/phyto/CoccoData/data/species_list.csv")
species_list = d['species']

m = library('/home/phyto/CoccoData/data/classification/phases.yml',
            '/home/phyto/CoccoData/data/classification/family.yml',
            "/home/phyto/CoccoData/data/sizes/",
            "/home/phyto/CoccoData/data/pic/",
            species_list)

#create a list of HOL species for future use:
HOL_list = m.return_HOL_list()

#return the named tuple so we can use it:
ntpl = m.return_ntpl()


lm_library_mean =  pd.Series([t.mean for t in ntpl  if t.id == "young2024" or t.id == "viliot2021"])
lm_library_spp =  pd.Series([t.species for t in ntpl  if t.id == "young2024" or t.id == "viliot2021"])
lm = pd.DataFrame({'species': lm_library_spp, 'mean': lm_library_mean})
lm = lm.groupby(['species']).mean()
lm['method'] = 'lm'

sem_library_mean =  [t.mean for t in ntpl  if t.id == "sheward2024"]
sem_library_spp =  [t.species for t in ntpl  if t.id == "sheward2024"]
sem = pd.DataFrame({'species': sem_library_spp, 'mean': sem_library_mean})
sem = sem.groupby(['species']).mean()
sem['method'] = 'sem'

d = pd.concat([sem, lm])
d.reset_index(inplace=True)

d = d.pivot(index='species', columns='method', values='mean')

fig, ax = plt.subplots()

d = d.dropna()
no_drop = len(d)
x = np.log(d['lm'])
y = np.log(d['sem'])
sns.scatterplot(x=x, y=y, ax=ax, color="firebrick", s=100)


d['ratio'] = d['sem']/d['lm']
d.reset_index(inplace=True)
d_copy = d.copy()

d1 = d_copy[d_copy['ratio']<0.2]
d2 = d_copy[d_copy['ratio']>5]
d_dropped = pd.concat([d1, d2])

d = d[d['ratio']>0.2]
d = d[d['ratio']<5]

yes_drop = len(d)




#sns.scatterplot(x=d['lm'], y=d['sem'])
#plt.show()

x = np.log(d['lm'])
y = np.log(d['sem'])
abline_plot(slope=1, intercept=0, ax=ax, color='black', linestyle='dashed')
sns.scatterplot(x=x, y=y, ax=ax, s=100)


slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("slope :" + str(slope))
print("intercept: " + str(intercept))
print("r_value: " + str(r_value))

d_dropped.reset_index(inplace=True)
# for i in range(0, len(d_dropped)):
#     ax.text(x.iloc[i], y.iloc[i], d_dropped['species'].iloc[i], 
#         bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})

ax.set_xlim(3, np.nanmax([x, y]))
ax.set_ylim(3, np.nanmax([x, y]))


fig.suptitle('Cell size estimates SEM vs LM', size=20,  weight='bold')

plt.show()

print("fin")
import numpy as np 


#G. flabellatus estimate:

pic_per_lith = [9.65, 11.62, 4.25] #Guerreiro2021
observations_per_estimate = [213, 106, 108] #Guerreiro2021
number_of_liths = [40, 86, 56, 18, 120]   #yang and wei 2003

p=observations_per_estimate/np.sum(observations_per_estimate) 

pic_sphere = np.random.choice(pic_per_lith, 10000, p=p)*np.random.choice(number_of_liths, 10000)

mean_pic_per_sphere =  np.round(np.mean(pic_sphere), -1)
sd_pic_per_sphere = np.round(np.std(pic_sphere), -1)

#round to two significant figures:



print("fin")



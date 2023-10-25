import numpy as np 


#Pontosphaera estimate:

pic_per_lith = [77.97, 133.07, 237.68] #Guerreiro2021
observations_per_estimate = [13, 5, 2] #Guerreiro2021
number_of_liths = [45]   #yang and wei 2003

p=observations_per_estimate/np.sum(observations_per_estimate) 

pic_sphere = np.random.choice(pic_per_lith, 10000, p=p)*np.random.choice(number_of_liths, 10000)

mean_pic_per_sphere =  np.round(np.mean(pic_sphere), -2)
sd_pic_per_sphere = np.round(np.std(pic_sphere), -2)

#round to two significant figures:

shape_factor = 0.05

print("fin")



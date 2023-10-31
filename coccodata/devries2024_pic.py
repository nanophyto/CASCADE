import numpy as np 


def estimate_pontosphaera():
    pic_per_lith = [77.97, 133.07, 237.68] #Guerreiro2021
    observations_per_estimate = [13, 5, 2] #Guerreiro2021
    number_of_liths = [45]   #yang and wei 2003

    p=observations_per_estimate/np.sum(observations_per_estimate) 

    pic_sphere = np.random.choice(pic_per_lith, 10000, p=p)*np.random.choice(number_of_liths, 10000)

    mean_pic_per_sphere =  np.round(np.mean(pic_sphere), -2)*0.12
    sd_pic_per_sphere = np.round(np.std(pic_sphere), -2)*0.12

    d = {'species':['undefined Pontosphaera'], 
        'mean': [mean_pic_per_sphere],
        'std': [sd_pic_per_sphere],
        'ref': ['guerreiro2021 and yang2003']}

    d = pd.DataFrame(d)

    return(d)


def estimate_flabellatus():

    pic_per_lith = [9.65, 11.62, 4.25] #Guerreiro2021
    observations_per_estimate = [213, 106, 108] #Guerreiro2021
    number_of_liths = [40, 86, 56, 18, 120]   #yang and wei 2003

    p=observations_per_estimate/np.sum(observations_per_estimate) 

    pic_sphere = np.random.choice(pic_per_lith, 10000, p=p)*np.random.choice(number_of_liths, 10000)

    mean_pic_per_sphere =  np.round(np.mean(pic_sphere), -1)*0.12
    sd_pic_per_sphere = np.round(np.std(pic_sphere), -1)*0.12


    d = {'species':['Gladiolithus flabellatus'], 
        'mean': [mean_pic_per_sphere],
        'std': [sd_pic_per_sphere],
        'ref': ['guerreiro2021 and yang2003']}

    d = pd.DataFrame(d)

d_flabellatus = estimate_flabellatus()
d_pontosphaera = estimate_pontosphaera()

d = pd.concat([d_flabellatus, d_pontosphaera])

d.to_csv("/home/phyto/CoccoData/pic/devries2024_pic.csv", index=False)


print("fin")



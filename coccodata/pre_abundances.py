

import pandas as pd
import numpy as np
from yaml import load, Loader

def preprocess_data():

    def clean_hagino200(): 
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Hagino2006.csv", na_values = "-")

        d['Latitude'] = d['Latitude'].str.rstrip('’S')
        d['Latitude'] = d['Latitude'].str.rstrip('’N')
        d['Latitude'] = d['Latitude'].str.replace('˚', '.')


        d['Longitude'] = d['Longitude'].str.rstrip('’E')
        d['Longitude'] = d['Longitude'].str.rstrip('’W')
        d['Longitude'] = d['Longitude'].str.replace('˚', '.')

        date = pd.to_datetime(d['Date'], format="%d/%m/%Y")

        d = d.groupby(by=d.columns, axis=1).sum()

        d['Reference'] = 'Hagino2006'
        d['Date'] = date
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns='Date', inplace=True)

        d = d.replace('--',np.NaN)

        d['Method'] = "SEM"
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.fillna(0, inplace=True)

        d = d.reset_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.to_csv("/home/phyto/CoccoData/data/abundances/Hagino2006.csv", index=False)



    def clean_devries2021(): 


        with open('/home/phyto/CoccoData/data/classification/pangaea.yml', 'r') as f:
            new_names = load(f, Loader=Loader)


        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/deVries-etal_2020.csv")

        d['Date/Time'] = pd.to_datetime(d['Date/Time'])

        d = d[d.Reference != 'Takahashi & Okada (2000)'] #wrong data

        d.rename(columns = {'Depth water [m]':'Depth'}, inplace = True)

        #remove samples below 200m (9 samples)
        d = d[d.Depth <= 201] 


        #Saavedra-Pellitero et al. (2014) samples:
        #weirdly lablled lon

        d = d.replace(-182.73, 2.73)
        d = d.replace(-182.85, 2.85)
        d = d.replace(-185.47, 5.47)


        d.rename(columns=new_names, inplace=True)
        d = d.groupby(by=d.columns, axis=1).sum()

        d.loc[d["Reference"] == "Saavedra-Pellitero et al. (2014)", "Gephyrocapsa ericsonii"] = 0

        d['Day'] = pd.DatetimeIndex(d['Date/Time']).day
        d['Month'] = pd.DatetimeIndex(d['Date/Time']).month
        d['Year'] = pd.DatetimeIndex(d['Date/Time']).year

        
        d.drop(columns=['Date/Time',"Sample method"], inplace=True)

        d.rename(columns = {'Depth water [m]':'Depth'}, inplace = True)

        d['Method'] = "SEM"


        with open('/home/phyto/CoccoData/data/unprocessed/references/refs_devries.yml', 'r') as f:
            short_refs = load(f, Loader=Loader)

        d=d.replace({"Reference": short_refs})

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.reset_index(inplace=True)

        for i in range(len(d['Reference'].unique())):
            reference_to_export = d['Reference'].unique()[i]
            d[d.Reference==reference_to_export].to_csv("/home/phyto/CoccoData/data/abundances/" + reference_to_export + ".csv", index=False)



    def clean_okada1973(): 


        d1 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/Okada1973.csv")
        d2 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/Okada1973B.csv")
        d3 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse2.csv")
        d4 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse4.csv")
        d5 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse4_A.csv")
        d6 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse5.csv")
        d7 = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/okada/traverse13.csv")

        d = pd.concat([d1, d2, d3, d4, d5, d6, d7])

        d['Method'] = "SEM"
        d['Reference'] = "Okada1973"
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.drop(['Temperature', 'station', 'Sample', 'Date'], axis = 1)

        def cell_per_litre(d):

            #find the number of cells counted:
            total = d['Total'].values

            d = d.drop(['Total'], axis = 1)
            counts= d.sum(axis=1).values

            #loop over each row and estimate the cells per liter by dividing by the number of counts and multiplying by the total:
            
            for x in range(0, d.shape[0]): #for each row:
                for i in range(0, d.shape[1]): #for each column:
                    d.iloc[x, i] = np.round((d.iloc[x, i]/counts[x]) * total[x])
            return(d)

        d = cell_per_litre(d)
        d = d.reset_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d.to_csv("/home/phyto/CoccoData/data/abundances/Okada1973.csv", index=False)



    def clean_cortes2001():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/cortes2001.csv", na_values = "ND")
        
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date'], inplace=True)
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Cortes2001.csv", index=False)


    def clean_takahashi2000():

        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Takahashi_raw.csv")
        d['Reference'] = "Takahashi2000"
        d['Method'] = "SEM"
        d['Date/Time'] = pd.to_datetime(d['Date/Time'])
        d['Day'] = pd.DatetimeIndex(d['Date/Time']).day
        d['Month'] = pd.DatetimeIndex(d['Date/Time']).month
        d['Year'] = pd.DatetimeIndex(d['Date/Time']).year

        counts =  d['counts'].values
        total = d['total'].values

        d = d.drop(['counts', 'total', 'Date/Time'], axis = 1)
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])

        for x in range(0, d.shape[0]): #for each row:
            for i in range(0, d.shape[1]): #for each column:
                d.iloc[x, i] = np.round((d.iloc[x, i]/counts[x]) * total[x])

        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Takahashi2010.csv", index=False)


    def clean_sal2013():

        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/sal2013.csv")
        d['Reference'] = "Sal2013"
        d.rename(columns = {'Lat':'Latitude'}, inplace = True)
        d.rename(columns = {'Lon':'Longitude'}, inplace = True)
        d['Method'] = "LM"

        d['Date'] = pd.to_datetime(d['Date'])
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d.drop(columns=['Date'], inplace=True)

        d = d[d['Latitude'].notna()]
        d = d[d['Longitude'].notna()]
        d = d[d['Depth'].notna()]
        d = d[d['Month'].notna()]
        d = d[d['Year'].notna()]
        d = d[d['Day'].notna()]

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Sal2013.csv", index=False)


    def clean_obrien2013():


        with open('/home/phyto/CoccoData/data/unprocessed/references/references_obrien.yml', 'r') as f:
            references = load(f, Loader=Loader)

        #reverse key-binding to: long-short
        references = {v: k for k, v in references.items()}


        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/obrien2013.csv")

        d = d.replace(["Utermohl", 
                "phase-contrast microscope",
                "Sedimentation and Inverted Microscope",
                    "Inverted Microscope",
                    "Light Microscope (x400)",
                    "epifluorescence microscopy",
                    "inverted microscope",
                    "light microscopy",
                "Inverted light and epifluorescent microscopy (Utermohl)",
                'Inverted microscope (Utermohl)',
                'Inverted microscope',
                ], "LM")

        print(d['Instrument/Method'].unique())

        d = d[d["Instrument/Method"].isin(["LM"])]

        d.rename(columns = {'Method References':'Reference'}, inplace = True)
        d.rename(columns = {'Depth (m)':'Depth'}, inplace = True)
        d.rename(columns = {'Longitutude':'Longitude'}, inplace = True)
        d.rename(columns = {'Instrument/Method':'Method'}, inplace = True)


        d = d.pivot_table(index=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'], columns=['Species/lowest classification'], values=['Cells/l'])
        d.columns = d.columns.droplevel(0)
        d = d.reset_index().rename_axis(None, axis=1)

        #d.reset_index(inplace=True)
        d=d.replace({"Reference": references})
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/OBrien2013.csv", index=False)



    def clean_estrada2016():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Malaspina2010.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day

        d.drop(columns=['Date', 'Station'], inplace=True)

        d['Method'] = "LM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Estrada2016.csv", index=False)



    def clean_baumann2000():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Baumann2000.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).year

        d.drop(columns=['Date'], inplace=True)
        d['Reference'] = "Baumann2000"
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Baumann2000.csv", index=False)


    def clean_kleijne1984():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Kleijne1984.csv")

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d.drop(columns=['Date'], inplace=True)
        d['Method'] = "SEM"

        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'])

        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Kleijne1984.csv", index=False)


    def clean_keuter2023():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Keuter2023.csv")
        d['Method'] = "SEM"
        d['Reference'] = "Keuter2023"

        d = d.set_index(['Latitude', 'Longitude', 'Depth','Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000
        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Keuter2023.csv", index=False)

    
    def clean_keuter2022():
        d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/abundances/Keuter2022.csv")
        d['Method'] = "SEM"
        d['Reference'] = "Keuter2022"

        d['Month'] = pd.DatetimeIndex(d['Date']).month
        d['Year'] = pd.DatetimeIndex(d['Date']).year
        d['Day'] = pd.DatetimeIndex(d['Date']).day
        d.drop(columns=['Date'], inplace=True)

        d = d.set_index(['Latitude', 'Longitude', 'Depth','Day', 'Month', 'Year', 'Reference', 'Method'])
        d = d*1000

        d = d.reset_index()
        d.to_csv("/home/phyto/CoccoData/data/abundances/Keuter2022.csv", index=False)
    
    clean_estrada2016()
    clean_hagino200()
    clean_devries2021()
    clean_okada1973()
    clean_takahashi2000()
    clean_cortes2001()
    clean_baumann2000()
    clean_sal2013()
    clean_obrien2013()
    clean_kleijne1984()
    clean_keuter2023()
    clean_keuter2022()

d = preprocess_data()


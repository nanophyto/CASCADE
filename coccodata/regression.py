import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt


class regression_simulation():
    """
    A function to simulate POC estimates based on regression.

    Parameters
    ----------
        X_train : float64
            Cell volume training data (um^3)
        X_predict : float64
            Cell volume used to estimate POC (um^3)
        Y_train : float64
            Carbon content training data (pg POC)
        n: int
            Number of random samples to simulate

    Returns
    -------
    simulated POC content (pg C) in array of length n

    Example
    -------
    d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
    Y_train = d['pg poc']
    X_train = d['volume']
    X_predict = [100, 99, 98, 101, 102, 100]*1000
    m = regression_simulation(X_train, X_predict, Y_train)
    test1 = m.simulate_data()
    sns.histplot(x=test1)
    plt.show()
    
    """
    def __init__(self, X_train, X_predict, Y_train):

        self.model = sm.GLM(Y_train, X_train, family=sm.families.Gamma(sm.families.links.identity())) # sm.OLS(Y, X)
        self.results = self.model.fit()
        self.X_predict = X_predict
        self.X_train = X_train
        self.Y_train = Y_train
    

    def regression_simulation_sample(self, X_predict_sample):
        
        gen = self.model.get_distribution(params = self.results.params, scale = self.results.scale, exog=X_predict_sample)

        simulated_data = gen.rvs(len(self.X_train))

        # simulated_data = simulated_data[simulated_data > 0] not needed for Gamma
        simulated_sample = np.random.choice(simulated_data, 1)

        return(simulated_sample)

    def simulate_data(self):

        simulated_data = []

        for i in range(len(self.X_predict)):
            simulated_data.extend(self.regression_simulation_sample(self.X_predict[i]))

        return(simulated_data)

# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
# Y_train = d['pg poc']
# X_train = d['volume']
# X_predict = [100, 99, 98, 101, 102, 100]*1000
# m = regression_simulation(X_train, X_predict, Y_train)
# test1 = m.simulate_data()
# sns.histplot(x=test1)
# plt.show()

print("fin")
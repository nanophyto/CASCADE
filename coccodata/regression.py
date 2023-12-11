import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt


def regression_simulation(X_train, X_predict, Y_train, n):
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
    Y = d['pg poc']
    X = d['volume']
    test1 = regression_simulation(X, 100, Y, 100)
    sns.histplot(x=test1)
    plt.show()
    
    """
    def regression_simulation_sample(X_train, X_predict, Y_train):
            
        model = sm.GLM(Y_train, X_train, family=sm.families.Gamma(sm.families.links.identity())) # sm.OLS(Y, X)
        results = model.fit()

        gen = model.get_distribution(params = results.params, scale = results.scale, exog=100)

        simulated_data = gen.rvs(len(X_train))
        # simulated_data = simulated_data[simulated_data > 0] not needed for Gamma
        simulated_sample = np.random.choice(simulated_data, 1)
        return(simulated_sample)


    simulated_data = []
    for i in range(n):
        simulated_data.extend(regression_simulation_sample(X_train, X_predict, Y_train,))

    return(simulated_data)


print("fin")
import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.graphics.api import abline_plot
from sklearn.metrics import r2_score
from sklearn.preprocessing import OneHotEncoder


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

        #Y_train = np.log10(Y_train)
        #X_train = np.log10(X_train)

        self.model = sm.GLM(Y_train, sm.add_constant(X_train), family=sm.families.Gamma(sm.families.links.identity())) # sm.OLS(Y, X)
        #self.model = sm.GLM(Y_train, X_train, family=sm.families.Gaussian()) # sm.OLS(Y, X)
        print(X_train)


        self.results = self.model.fit()
        self.X_predict = X_predict
        self.X_train = X_train
        self.Y_train = Y_train
    
    def regression_simulation_sample(self, X_predict_sample):
    
        gen = self.model.get_distribution(params = self.results.params, scale = self.results.scale, exog=X_predict_sample)
        simulated_data = gen.rvs()
        # simulated_data = simulated_data[simulated_data > 0] not needed for Gamma
        simulated_sample = np.random.choice(simulated_data, 1)

        return(simulated_sample)

    def simulate_data(self):

        simulated_data = []

        for i in range(len(self.X_predict)):
            simulated_data.extend(self.regression_simulation_sample(self.X_predict[i]))

        return(simulated_data)


    def return_performance(self):
        print(self.results.summary())
        RMSE = sm.tools.eval_measures.rmse(self.Y_train, self.results.mu, axis=0)
        print("RMSE: " + str(RMSE))
        MAE = sm.tools.eval_measures.meanabs(self.Y_train, self.results.mu, axis=0) 
        print("MAE: " + str(MAE))

        rRMSE = RMSE/np.mean(self.Y_train)
        rMAE = MAE/np.mean(self.Y_train)
        print("mean:" + str(np.mean(self.Y_train)))
        print("rRMSE: " + str(rRMSE))
        print("rMAE: " + str(rMAE))
        aic = self.results.aic
        print("AIC: " + str(aic))

        bias = sm.tools.eval_measures.meanabs(self.Y_train, self.results.mu, axis=0) 
        print("bias: " + str(bias))

    def plot_fit(self):
        y = self.Y_train
        yhat = self.results.mu
        # log transform

        y = np.log10(y)
        yhat = np.log10(yhat)
        fig, ax = plt.subplots()
        ax.scatter(yhat, y)
        line_fit = sm.OLS(y, sm.add_constant(yhat, prepend=True)).fit()
        abline_plot(model_results=line_fit, ax=ax)

        abline_plot(slope=1, intercept=0, ax=ax, color='black', linestyle='dashed')


        ax.set_title('Model Fit Plot')
        ax.set_ylabel('Observed values (pg C, log10)')
        ax.set_xlabel('Fitted values (pg C, log10)')
        #ax.set_xlim(0, np.nanmax([y, yhat]))
        #ax.set_ylim(0, np.nanmax([y, yhat]))

        plt.show()


# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
# Y_train = d['pg poc']
# X_train = d['volume']

# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/rosie_size_pic.csv")
# d = d.dropna()
# d = d[d['PIC pg C'] >0]
# d = d[d['Volume'] >0]
# Y_train = d['PIC pg C']

# d = pd.concat([pd.get_dummies(d, columns=['Phase'], dtype=float), d['Volume']], axis=1)
# X_train = d[['Volume', 'Phase_HET', 'Phase_HOL']].astype(float)
# X_train = d['Volume']

# X_predict = [100, 99, 98, 101, 102, 100]*1000
# m = regression_simulation(X_train, X_predict, Y_train)
# # test1 = m.simulate_data()
# # sns.histplot(x=test1)
# # plt.show()
# m.return_performance()
# m.plot_fit()
# print("fin")

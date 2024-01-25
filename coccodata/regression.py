import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.graphics.api import abline_plot
from sklearn.metrics import r2_score
from sklearn.preprocessing import OneHotEncoder
from scipy import stats
plt.rcParams.update({'font.size': 14})
from functions import bootstrap


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

        X_train = bootstrap(X_train)
        Y_train = bootstrap(Y_train)
        self.model = sm.GLM(Y_train,  X_train, family=sm.families.Gamma(sm.families.links.identity())) # sm.OLS(Y, X)
        #self.model = sm.GLM(Y_train, X_train, family=sm.families.Gaussian()) # sm.OLS(Y, X)
        #print(X_train)

        self.results = self.model.fit()
        self.X_predict = X_predict
        self.X_train = X_train
        self.Y_train = Y_train
    
    def regression_simulation_sample(self, X_predict_sample):
    
        gen = self.model.get_distribution(params = self.results.params, scale = self.results.scale, exog=X_predict_sample)
        simulated_data = np.atleast_1d(gen.rvs())
        #simulated_data = simulated_data[simulated_data > 0] #not needed for Gamma
        #simulated_sample = np.random.choice(simulated_data, 1)
        #simulated_data = bootstrap(simulated_data)
        return(simulated_data)

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

    def plot_fit_POC(self):
        y = self.Y_train
        x = self.X_train
        yhat = self.results.mu
        # log transform
        y = np.log10(y)
        yhat = np.log10(yhat)
        x = np.log10(x)


        fig, axs = plt.subplots(1, 2)

        axs[0].scatter(x, y)
        sns.regplot(x=x, y=y, order=1, ax=axs[0], ci=None, scatter=False)

        axs[1].scatter(yhat, y)

        line_fit = sm.OLS(y, sm.add_constant(yhat, prepend=True)).fit()
        abline_plot(model_results=line_fit, ax=axs[1])
        abline_plot(slope=1, intercept=0, ax=axs[1], color='black', linestyle='dashed')


        axs[0].set_title('Allometric scaling')
        axs[0].set_ylabel('Carbon content (pg C, log10)')
        axs[0].set_xlabel('Cell size (um3, log10)')

        axs[1].set_title('Observed vs Predicted')
        axs[1].set_ylabel('Observed values (pg C, log10)')
        axs[1].set_xlabel('Predicted values (pg C, log10)')

        axs[1].set_xlim(0, np.nanmax([y, yhat]))

        axs[0].text(0.05, 0.95, "A)", transform=axs[0].transAxes,
            fontsize=16, fontweight='bold', va='top')
        axs[1].text(0.05, 0.95, "B)", transform=axs[1].transAxes,
            fontsize=16, fontweight='bold', va='top')

        fig.suptitle('Allometric GLM for POC', size=24,  weight='bold')

        plt.show()

        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        print("slope :" + str(slope))
        print("intercept: " + str(intercept))
        print("r_value: " + str(r_value))



    def plot_fit_PIC(self):
        
        y = self.Y_train
        x = self.X_train['Volume']
        yhat = self.results.mu
        # log transform
        y = np.log10(y)
        yhat = np.log10(yhat)
        x = np.log10(x)


        y_HOL = y[self.X_train['Phase_HOL']==1]
        yhat_HOL = yhat[self.X_train['Phase_HOL']==1]

        y_HET = y[self.X_train['Phase_HET']==1]
        yhat_HET = yhat[self.X_train['Phase_HET']==1]


        x_HET = x[self.X_train['Phase_HET']==1]
        x_HOL = x[self.X_train['Phase_HOL']==1]

        HOL_line_fit = sm.OLS(y_HOL, sm.add_constant(yhat_HOL, prepend=True)).fit()
        HET_line_fit = sm.OLS(y_HET, sm.add_constant(yhat_HET, prepend=True)).fit()



        fig, axs = plt.subplots(2,2)

        #scaling HET
        axs[0,0].scatter(x_HET, y_HET, color='firebrick')
        sns.regplot(x=x_HET, y=y_HET, order=1, ax=axs[0,0], color='darkred', ci=None, scatter=False)
        axs[0,0].set_title('Scaling (diploid)')
        axs[0,0].set_ylabel('C content (pg C)')
        axs[0,0].set_xlabel('Cell size (um3)')


        #Fit HET
        axs[0,1].scatter(yhat_HET, y_HET, color='firebrick')
        abline_plot(model_results=HET_line_fit, ax=axs[0,1], color='darkred')
        abline_plot(slope=1, intercept=0, ax=axs[0,1], color='black', linestyle='dashed')
        axs[0,1].set_title('Fit (diploid)')
        axs[0,1].set_ylabel('Observed (pg C)')
        axs[0,1].set_xlabel('Predicted (pg C)')
        axs[0,1].set_xlim(0, np.nanmax([y, yhat]))
        axs[0,1].set_ylim(0, np.nanmax([y, yhat]))



        # scaling HOL
        axs[1,0].scatter(x_HOL, y_HOL, color='steelblue')
        sns.regplot(x=x_HOL, y=y_HOL, order=1, ax=axs[1,0], color='darkblue', ci=None, scatter=False)
        axs[1,0].set_title('Scaling (haploid)')
        axs[1,0].set_ylabel('C content (pg C)')
        axs[1,0].set_xlabel('Cell size (um3)')


        #Fit HOL
        axs[1,1].scatter(yhat_HOL, y_HOL, color='steelblue')
        abline_plot(model_results=HOL_line_fit, ax=axs[1,1], color='darkblue')
        abline_plot(slope=1, intercept=0, ax=axs[1,1], color='black', linestyle='dashed')
        axs[1,1].set_title('Fit (haploid)')
        axs[1,1].set_ylabel('Observed (pg C)')
        axs[1,1].set_xlabel('Predicted (pg C)')
        axs[1,1].set_xlim(0, np.nanmax([y, yhat]))
        axs[1,1].set_ylim(0, np.nanmax([y, yhat]))

        # axs[1,0].set_aspect('equal', adjustable='box')
        # axs[1,1].set_aspect('equal', adjustable='box')
        # axs[0,0].set_aspect('equal', adjustable='box')
        # axs[0,1].set_aspect('equal', adjustable='box')

        axs[0,0].text(0.05, 0.95, "A)", transform=axs[0,0].transAxes,
            fontsize=16, fontweight='bold', va='top')
        axs[0,1].text(0.05, 0.95, "B)", transform=axs[0,1].transAxes,
            fontsize=16, fontweight='bold', va='top')
        axs[1,0].text(0.05, 0.95, "C)", transform=axs[1,0].transAxes,
            fontsize=16, fontweight='bold', va='top')
        axs[1,1].text(0.05, 0.95, "D)", transform=axs[1,1].transAxes,
            fontsize=16, fontweight='bold', va='top')

        fig.suptitle('Allometric GLM for PIC', size=24,  weight='bold')


        plt.tight_layout()
        plt.show()

# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
# Y_train = d['pg poc']
# X_train = d['volume']
# X_predict = np.random.normal(8, 2, 1000)
# m = regression_simulation(X_train, X_predict, Y_train)
# test1 = m.simulate_data()
# sns.histplot(x=test1)
# plt.show()
# m.return_performance()

# m.plot_fit_POC()



# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/rosie_size_pic.csv")
# d = d.dropna()
# d = d[d['PIC pg C'] >0]
# d = d[d['Volume'] >0]
# Y_train = d['PIC pg C']

# d = pd.get_dummies(d, columns=['Phase'], dtype=float)
# X_train = d[['Volume', 'Phase_HET', 'Phase_HOL']].astype(float)
# #X_train = d['Volume']

# X_predict = [100, 99, 98, 101, 102, 100]*1000
# m = regression_simulation(X_train, X_predict, Y_train)
# # test1 = m.simulate_data()
# # sns.histplot(x=test1)
# # plt.show()
# m.return_performance()
# m.plot_fit_PIC()


print("fin")

import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
from joblib import parallel_backend, Parallel, delayed
from statsmodels.graphics.api import abline_plot
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import OneHotEncoder
from scipy import stats

plt.rcParams.update({"font.size": 14})
from functions import bootstrap


class regression_simulation:
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
        # won't this break ordering between Y and x?
        # see @ljwolf's implementation later, in the
        # regression coefficient example
        allometric_model = sm.GLM(
            Y_train,
            np.log(X_train),
            family=sm.families.Gamma(link=sm.families.links.Log()),
        )
        # self.model = sm.GLM(Y_train, X_train, family=sm.families.Gaussian()) # sm.OLS(Y, X)
        # print(X_train)

        self._fit_function = lambda y, x: (
            sm.GLM(
                y, np.log(x), family=sm.families.Gamma(link=sm.families.links.Log())
            ).fit()
        )
        self.results = allometric_model.fit()
        self.X_predict = X_predict
        self.X_train = X_train
        self.Y_train = Y_train

    def simulate_regression_params(self, boot=False, n_replicates=1000):
        """
        Simulate a slope (and potentially intercept) used
        for the construction of prediction intervals and
        replicate distributions.

        Parameters
        ----------
        boot : bool (Default: False)
            whether or not to use bootstrapping on the
            fitted data to simulate regression parameters.
            If not (default), we use the MLE asymptotic
            behavior-parameters are normally distributed
            with a covariance matrix from the variance-covariance
            matrix in the model.
        n_replicates : int (Default: 1000)
            how many replications of the regression parameters
            to draw.

        Returns
        -------
        numpy.ndarray of shape (n_replicates, n_params), where
        n_params is 1 if there is only a slope fit in the
        initial model, and >2 if there is both a slope
        and an intercept or multiple slope parameters.
        """
        if not boot:
            # this leverages the fact that the regression
            # MLE is asymptotically gaussian. That is,
            # the maximum likelihood estimates of
            # our parameters (slope, intercept)
            # will be asymptotically gaussian with a mean
            # of their estimated values and a variance-covariance
            # matrix according to (x'x)^{-1}sigma^2
            gradient = self.results.params
            gradient_uncertainty = self.results.cov_params()
            if len(gradient) == 0:
                # this is one slope with one variance
                return np.random.normal(
                    gradient, gradient_uncertainty, size=(n_replicates, 1)
                )
            else:
                # this is possibly an intercept and a slope,
                # with a variance-covariance matrix, showing
                # how the uncertainty is correlated across
                # parameters
                return np.random.multivariate_normal(
                    gradient, gradient_uncertainty, size=n_replicates
                )
        else:
            # If we don't use the mle asymptotics,
            # we bootstrap the input data, re-estimate
            # the slope, and use that as a "replicate"
            # of the possible slope.
            df = pd.DataFrame(
                np.column_stack(
                    (self.results.model.data.endog, self.results.model.data.exog)
                ),
                columns=["y", *self.results.model.data.cov_names],
            )

            def replicate(_, df=df):
                """
                generate one bootstrapped replication
                using the inputted data
                """
                data = df.sample(frac=1, replace=True)
                return self._fit_function(data.y, np.exp(data.drop("y", axis=1))).params

            engine = Parallel()
            promise = delayed(replicate)
            return np.row_stack(engine(promise(i) for i in range(n_replicates)))

    def simulate_data(self):
        simulated_data = []

        for i in range(len(self.X_predict)):
            simulated_data.extend(self.regression_simulation_sample(self.X_predict[i]))

        return simulated_data

    def simulate_predictions(
        self,
        X_predict=None,
        params=None,
        conditional=False,
        n_replicates=1,
        boot=None,
    ):
        """
        Simulate from a predictive distribution.

        Parameters
        ----------
        X_predict : numpy.ndarray (default: None)
            values corresponding to the x values where
            we want to generate a distribution of predictions.
            By default, this uses the locations provided
            to the class upon initialisation, in self.X_predict
        params : numpy.ndarray (default: None)
            parameter values to use when simulating
            from the predictive distribution. If None,
            we use the fitted values from the simulation
            object. Must be None if the simulations are
            unconditional.
        conditional : bool (default: False)
            whether to condition the simulations on a "known"
            slope, removing the variation due to estimation
            error in the relationship between size and carbon
        n_replicates : int (default: 1)
            how many replicates to create for the locations
            in `X_predict`
        boot : bool (default: False)
            whether or not to use bootstrapping to simulate
            regression parameters. Ignored if using conditional
            simulation.

        Returns
        -------
        numpy.ndarray of size (X_predict, n_replicates) for
        prediction values at each inputted X_predict.
        """
        if X_predict is None:
            X_predict = np.atleast_2d(self.X_predict).reshape(-1, 1)
        if conditional:
            assert boot is False, "cannot bootstrap conditional predictions"
            if params is None:
                params = self.results.params
            gen = self.results.model.get_distribution(
                params=params,
                scale=self.results.scale,
                exog=X_predict,
            )
            simulated_data = np.column_stack(
                [np.atleast_1d(gen.rvs()) for _ in range(n_replicates)]
            )
        else:
            if params is not None:
                warnings.warn(
                    "location parameters (`params`)  provided, but unconditional distribution was requested. Ignoring the provided location parameters...",
                    stacklevel=2,
                )
            params = self.simulate_regression_params(
                n_replicates=n_replicates, boot=boot
            )
            mus = X_predict @ params
            simulated_data = np.column_stack(
                [self.results.model.family.predict(mu_vec) for mu_vec in mus.T]
            )
        return simulated_data

    def return_performance(self):
        print(self.results.summary())
        RMSE = sm.tools.eval_measures.rmse(self.Y_train, self.results.predict(), axis=0)
        print("RMSE: " + str(RMSE))
        MAE = sm.tools.eval_measures.meanabs(
            self.Y_train, self.results.predict(), axis=0
        )
        print("MAE: " + str(MAE))

        rRMSE = RMSE / np.mean(self.Y_train)
        rMAE = MAE / np.mean(self.Y_train)
        print("mean:" + str(np.mean(self.Y_train)))
        print("rRMSE: " + str(rRMSE))
        print("rMAE: " + str(rMAE))
        aic = self.results.aic
        print("AIC: " + str(aic))

        bias = sm.tools.eval_measures.meanabs(
            self.Y_train, self.results.predict(), axis=0
        )
        print("bias: " + str(bias))

    def plot_fit_POC(self):
        y = self.Y_train
        x = self.X_train
        yhat = self.results.predict()
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
        abline_plot(slope=1, intercept=0, ax=axs[1], color="black", linestyle="dashed")

        axs[0].set_title("Allometric scaling")
        axs[0].set_ylabel("Carbon content (pg C, log10)")
        axs[0].set_xlabel("Cell size (um3, log10)")

        axs[1].set_title("Observed vs Predicted")
        axs[1].set_ylabel("Observed values (pg C, log10)")
        axs[1].set_xlabel("Predicted values (pg C, log10)")

        axs[1].set_xlim(0, np.nanmax([y, yhat]))

        axs[0].text(
            0.05,
            0.95,
            "A)",
            transform=axs[0].transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )
        axs[1].text(
            0.05,
            0.95,
            "B)",
            transform=axs[1].transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )

        fig.suptitle("Allometric GLM for POC", size=24, weight="bold")

        plt.show()

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        print("slope :" + str(slope))
        print("intercept: " + str(intercept))
        print("r_value: " + str(r_value))

    def plot_fit_PIC(self):
        y = self.Y_train
        x = self.X_train["Volume"]
        yhat = self.results.mu
        # log transform
        y = np.log10(y)
        yhat = np.log10(yhat)
        x = np.log10(x)

        y_HOL = y[self.X_train["Phase_HOL"] == 1]
        yhat_HOL = yhat[self.X_train["Phase_HOL"] == 1]

        y_HET = y[self.X_train["Phase_HET"] == 1]
        yhat_HET = yhat[self.X_train["Phase_HET"] == 1]

        x_HET = x[self.X_train["Phase_HET"] == 1]
        x_HOL = x[self.X_train["Phase_HOL"] == 1]

        HOL_line_fit = sm.OLS(y_HOL, sm.add_constant(yhat_HOL, prepend=True)).fit()
        HET_line_fit = sm.OLS(y_HET, sm.add_constant(yhat_HET, prepend=True)).fit()

        fig, axs = plt.subplots(2, 2)

        # scaling HET
        axs[0, 0].scatter(x_HET, y_HET, color="firebrick")
        sns.regplot(
            x=x_HET,
            y=y_HET,
            order=1,
            ax=axs[0, 0],
            color="darkred",
            ci=None,
            scatter=False,
        )
        axs[0, 0].set_title("Scaling (diploid)")
        axs[0, 0].set_ylabel("C content (pg C)")
        axs[0, 0].set_xlabel("Cell size (um3)")

        # Fit HET
        axs[0, 1].scatter(yhat_HET, y_HET, color="firebrick")
        abline_plot(model_results=HET_line_fit, ax=axs[0, 1], color="darkred")
        abline_plot(
            slope=1, intercept=0, ax=axs[0, 1], color="black", linestyle="dashed"
        )
        axs[0, 1].set_title("Fit (diploid)")
        axs[0, 1].set_ylabel("Observed (pg C)")
        axs[0, 1].set_xlabel("Predicted (pg C)")
        axs[0, 1].set_xlim(0, np.nanmax([y, yhat]))
        axs[0, 1].set_ylim(0, np.nanmax([y, yhat]))

        # scaling HOL
        axs[1, 0].scatter(x_HOL, y_HOL, color="steelblue")
        sns.regplot(
            x=x_HOL,
            y=y_HOL,
            order=1,
            ax=axs[1, 0],
            color="darkblue",
            ci=None,
            scatter=False,
        )
        axs[1, 0].set_title("Scaling (haploid)")
        axs[1, 0].set_ylabel("C content (pg C)")
        axs[1, 0].set_xlabel("Cell size (um3)")

        # Fit HOL
        axs[1, 1].scatter(yhat_HOL, y_HOL, color="steelblue")
        abline_plot(model_results=HOL_line_fit, ax=axs[1, 1], color="darkblue")
        abline_plot(
            slope=1, intercept=0, ax=axs[1, 1], color="black", linestyle="dashed"
        )
        axs[1, 1].set_title("Fit (haploid)")
        axs[1, 1].set_ylabel("Observed (pg C)")
        axs[1, 1].set_xlabel("Predicted (pg C)")
        axs[1, 1].set_xlim(0, np.nanmax([y, yhat]))
        axs[1, 1].set_ylim(0, np.nanmax([y, yhat]))

        # axs[1,0].set_aspect('equal', adjustable='box')
        # axs[1,1].set_aspect('equal', adjustable='box')
        # axs[0,0].set_aspect('equal', adjustable='box')
        # axs[0,1].set_aspect('equal', adjustable='box')

        axs[0, 0].text(
            0.05,
            0.95,
            "A)",
            transform=axs[0, 0].transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )
        axs[0, 1].text(
            0.05,
            0.95,
            "B)",
            transform=axs[0, 1].transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )
        axs[1, 0].text(
            0.05,
            0.95,
            "C)",
            transform=axs[1, 0].transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )
        axs[1, 1].text(
            0.05,
            0.95,
            "D)",
            transform=axs[1, 1].transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )

        fig.suptitle("Allometric GLM for PIC", size=24, weight="bold")

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    d = pd.read_csv("../data/unprocessed/poulton2024.csv")
    Y_train = d["pg poc"]
    X_train = d["volume"]
    X_predict = np.random.normal(8, 2, 1000)

    loggamma_raw_x = sm.GLM(
        Y_train, X_train, family=sm.families.Gamma(link=sm.families.links.Log())
    ).fit()
    loggamma = sm.GLM(
        Y_train,
        np.log(X_train),
        family=sm.families.Gamma(link=sm.families.links.Log()),
    ).fit()
    ident_gamma = sm.GLM(
        Y_train, X_train, family=sm.families.Gamma(link=sm.families.links.Identity())
    ).fit()
    loglinear = sm.OLS(np.log(Y_train), np.log(X_train)).fit()

    support = np.linspace(X_train.min(), X_train.max())

    preds_gamma_raw_x = loggamma_raw_x.predict(support)
    preds_log_gamma = loggamma.predict(np.log(support))
    preds_ident_gamma = ident_gamma.predict(support)
    preds_loglinear = np.exp(loglinear.predict(np.log(support)))

    f, ax = plt.subplots(1, 1)
    ax.scatter(d["volume"], d["pg poc"], color="k")
    for pred_type, preds in [
        (k, v) for k, v in locals().items() if k.startswith("preds_")
    ]:
        if pred_type.endswith("raw_x"):
            continue
        ax.plot(
            support, preds, label=pred_type.lstrip("preds_").replace("_", " ").title()
        )
    ax.legend()
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.draw()

    m = regression_simulation(X_train, X_predict, Y_train)

    normal_approx = m.simulate_regression_params()
    booted_slopes = m.simulate_regression_params(boot=True)
    plot_data = pd.concat(
        (
            pd.DataFrame(normal_approx, columns=["slope estimate"]).assign(
                style="mle asymptotics"
            ),
            pd.DataFrame(booted_slopes, columns=["slope estimate"]).assign(
                style="boostrapped"
            ),
        )
    )
    f, ax = plt.subplots(1, 1, figsize=(6, 4))
    sns.kdeplot(plot_data, ax=ax, hue="style", x="slope estimate")
    plt.draw()

    conditional = m.simulate_predictions()
    unconditional = m.simulate_predictions(conditional=False)
    unconditional_booted = m.simulate_predictions(boot=True, conditional=False)

    plot_data = pd.concat(
        (
            pd.DataFrame(conditional, columns=["conditional"]),
            pd.DataFrame(unconditional, columns=["unconditional"]),
            pd.DataFrame(unconditional_booted, columns=["unconditional_booted"]),
        ),
        axis=1,
    ).melt(var_name="simulation_type", value_name="replicate")
    ax = sns.kdeplot(data=plot_data, hue="simulation_type", x="replicate")
    ax.set_xscale("log")
    plt.show()

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

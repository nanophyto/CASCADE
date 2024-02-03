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
        self.X_train = X_train
        self.X_predict = X_predict
        self.Y_train = Y_train
        allometric_model = sm.GLM(
            self.Y_train,
            self.X_train,
            family=sm.families.Gamma(link=sm.families.links.Log()),
        )
        self.results = allometric_model.fit()

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
                Y_train_ = data["y"].values
                X_train_ = data[self.results.model.data.cov_names]
                return (
                    sm.GLM(
                        Y_train_,
                        X_train_,
                        family=sm.families.Gamma(link=sm.families.links.Log()),
                    )
                    .fit()
                    .params
                )

            engine = Parallel()
            promise = delayed(replicate)
            return np.row_stack(engine(promise(i) for i in range(n_replicates)))

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
            if isinstance(self.X_predict, pd.DataFrame):
                ncol_X_pred = len(self.X_predict.columns)
            else:
                ncol_X_pred = 1
            X_predict = np.atleast_2d(self.X_predict).reshape(-1, ncol_X_pred)

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
            #            lin_preds = X_predict @ params.T

            lin_preds = X_predict @ params.T

            simulated_data = np.column_stack(
                [self.results.model.family.fitted(lin_pred) for lin_pred in lin_preds.T]
            )
        return simulated_data  # [simulated_data>0]

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

    def plot_fit_POC(self, x, y, boot=False, n_replicates=1, figsize=(12, 8)):
        # y = self.Y_train
        # x = self.X_train

        m = regression_simulation(x, x, y)
        yhat = pd.DataFrame(
            {
                "yhat": m.simulate_predictions(
                    boot=boot, n_replicates=n_replicates
                ).flatten()
            }
        )

        d = pd.concat([x, y, yhat], axis=1)

        fig, axs = plt.subplots(1, 2, figsize=figsize)

        axs[0].scatter(x, y)
        sns.regplot(
            x=d["volume"], y=d["pg poc"], order=1, ax=axs[0], ci=None, scatter=False
        )

        sns.regplot(x=d["yhat"], y=d["pg poc"], order=1, ax=axs[1], ci=None)
        abline_plot(slope=1, intercept=0, ax=axs[1], color="black", linestyle="dashed")

        axs[0].set_title("Allometric scaling")
        axs[0].set_ylabel("Carbon content (pg C, log10)")
        axs[0].set_xlabel("Cell size (um3, log10)")

        axs[1].set_title("Observed vs Predicted")
        axs[1].set_ylabel("Observed values (pg C, log10)")
        axs[1].set_xlabel("Predicted values (pg C, log10)")

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

    def plot_fit_PIC(
        self, figsize=(8, 8), title=None, log_trans=True, boot=False, n_replicates=1
    ):
        y = self.Y_train
        x = self.X_train  # ["Volume"]
        m = regression_simulation(x, x, y)
        yhat = pd.DataFrame(
            {
                "yhat": m.simulate_predictions(
                    boot=boot, n_replicates=n_replicates
                ).flatten()
            }
        )

        d = pd.concat([x, y, yhat], axis=1)
        d.dropna(inplace=True)
        # print(d)

        fig, axs = plt.subplots(2, 2, figsize=figsize)

        sns.scatterplot(
            data=d[d["Phase_HET"] == 1], x="Volume", y="PIC pg C", ax=axs[0, 0]
        )

        sns.scatterplot(
            data=d[d["Phase_HET"] == 1], y="PIC pg C", x="yhat", ax=axs[0, 1]
        )

        sns.scatterplot(
            data=d[d["Phase_HOL"] == 1], x="Volume", y="PIC pg C", ax=axs[1, 0]
        )

        sns.scatterplot(
            data=d[d["Phase_HOL"] == 1], y="PIC pg C", x="yhat", ax=axs[1, 1]
        )
        axs[0, 1].plot(
            d[d["Phase_HET"] == 1]["PIC pg C"], d[d["Phase_HET"] == 1]["PIC pg C"]
        )
        axs[1, 1].plot(
            d[d["Phase_HOL"] == 1]["PIC pg C"], d[d["Phase_HOL"] == 1]["PIC pg C"]
        )

        if log_trans == True:
            axs[0, 0].set_yscale("log")
            axs[0, 0].set_xscale("log")
            axs[0, 1].set_yscale("log")
            axs[0, 1].set_xscale("log")
            axs[1, 0].set_yscale("log")
            axs[1, 0].set_xscale("log")
            axs[1, 1].set_yscale("log")
            axs[1, 1].set_xscale("log")

        axs[0, 0].set_title("Scaling (diploid)")
        axs[0, 0].set_ylabel("C content (pg C)")
        axs[0, 0].set_xlabel("Cell size (um3)")

        axs[0, 1].set_title("Fit (diploid)")
        axs[0, 1].set_ylabel("Observed (pg C)")
        axs[0, 1].set_xlabel("Predicted (pg C)")

        axs[1, 0].set_title("Scaling (haploid)")
        axs[1, 0].set_ylabel("C content (pg C)")
        axs[1, 0].set_xlabel("Cell size (um3)")

        axs[1, 1].set_title("Fit (haploid)")
        axs[1, 1].set_ylabel("Observed (pg C)")
        axs[1, 1].set_xlabel("Predicted (pg C)")

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

        # fig.suptitle("Allometric GLM for PIC", size=24, weight="bold")
        # # axs[1,1].set_xscale("log")
        # # axs[0,1].set_xscale("log")
        # # axs[1,1].set_yscale("log")
        # # axs[0,1].set_yscale("log")

        # # axs[1,0].set_xscale("log")
        # # axs[0,0].set_xscale("log")
        # # axs[1,0].set_yscale("log")
        # # axs[0,0].set_yscale("log")
        fig.suptitle(title)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    test_POC = True
    test_PIC = True
    plot_PIC = True

    if test_POC == True:
        print("testing POC")
        #    d = pd.read_csv("../data/unprocessed/poulton2024.csv")
        d = pd.read_csv("../data/unprocessed/poulton2024.csv")
        Y_train = d["pg poc"]
        X_train = np.log(d["volume"] + 1)

        def diameter_to_volume(d):
            v = (1 / 6) * np.pi * d**3
            return v

        v = diameter_to_volume(10)
        X_predict = np.random.normal(v, v * 0.2, 1000)
        X_predict = np.log(X_predict[X_predict > 0] + 1)
        m = regression_simulation(X_train, X_predict, Y_train)

        poc = m.simulate_predictions(boot=True, n_replicates=1)
        median = np.percentile(poc, 50)
        print("boot = True:")
        print(median)
        print("expected value ~70")
        m.plot_fit_POC(X_train, Y_train, boot=True, n_replicates=1)

        poc = m.simulate_predictions(boot=False, n_replicates=1)
        median = np.percentile(poc, 50)
        print("boot = False:")
        print(median)
        print("expected value ~70")
        m.plot_fit_POC(X_train, Y_train, boot=False, n_replicates=1)

    if test_PIC == True:
        print("testing PIC")

        # d = pd.read_csv("../data/unprocessed/rosie_size_pic.csv")
        d = pd.read_csv("../data/unprocessed/rosie_size_pic.csv")
        d = d.dropna()
        d = d[d["PIC pg C"] > 0]
        d = d[d["Volume"] > 0]
        Y_train = d["PIC pg C"]

        def diameter_to_volume(d):
            v = (1 / 6) * np.pi * d**3
            return v

        v = diameter_to_volume(10)

        d = pd.get_dummies(d, columns=["Phase"], dtype=float)
        X_train = d[["Volume", "Phase_HET", "Phase_HOL"]].astype(float)
        X_train["Volume"] = np.log(X_train.Volume + 1)

        X_predict = pd.DataFrame(
            {
                "Volume": np.random.normal(v, v * 0.2, 1000),
                "Phase_HET": [1] * 1000,
                "Phase_HOL": [0] * 1000,
            }
        )
        X_predict = X_predict[X_predict["Volume"] > 0]
        X_predict["Volume"] = np.log(X_predict.Volume + 1)

        m = regression_simulation(X_train, X_predict, Y_train)

        pic = m.simulate_predictions(boot=False, n_replicates=1)
        median = np.percentile(pic, 50)
        print("boot = False:")
        print(median)
        print("expected value ~20")
        m.plot_fit_PIC(title="conditional", boot=False, n_replicates=1)
        # m.return_performance()

        # m.return_performance()
        # if plot_PIC == True:
        #    sns.histplot(x=pic.flatten())
        # plt.show()

        m.plot_fit_PIC(title="boot", boot=True, n_replicates=1)
        pic = m.simulate_predictions(boot=True, n_replicates=1)
        median = np.percentile(pic, 50)
        print("boot = True:")
        print(median)
        print("expected value ~20")
        # m.return_performance()

    # d = pd.read_csv("../data/unprocessed/poulton2024.csv")
    # Y_train = d["pg poc"]
    # X_train = d["volume"]
    # X_predict = np.random.normal(8, 2, 1000)

    # loggamma_raw_x = sm.GLM(
    #     Y_train, X_train, family=sm.families.Gamma(link=sm.families.links.Log())
    # ).fit()
    # loggamma = sm.GLM(
    #     Y_train,
    #     np.log(X_train),
    #     family=sm.families.Gamma(link=sm.families.links.Log()),
    # ).fit()
    # ident_gamma = sm.GLM(
    #     Y_train, X_train, family=sm.families.Gamma(link=sm.families.links.Identity())
    # ).fit()
    # loglinear = sm.OLS(np.log(Y_train), np.log(X_train)).fit()

    # support = np.linspace(X_train.min(), X_train.max())

    # preds_gamma_raw_x = loggamma_raw_x.predict(support)
    # preds_log_gamma = loggamma.predict(np.log(support))
    # preds_ident_gamma = ident_gamma.predict(support)
    # preds_loglinear = np.exp(loglinear.predict(np.log(support)))

    # f, ax = plt.subplots(1, 1)
    # ax.scatter(d["volume"], d["pg poc"], color="k")
    # for pred_type, preds in [
    #     (k, v) for k, v in locals().items() if k.startswith("preds_")
    # ]:
    #     if pred_type.endswith("raw_x"):
    #         continue
    #     ax.plot(
    #         support, preds, label=pred_type.lstrip("preds_").replace("_", " ").title()
    #     )
    # ax.legend()
    # ax.set_yscale("log")
    # ax.set_xscale("log")
    # plt.draw()

    # m = regression_simulation(X_train, X_predict, Y_train)

    # normal_approx = m.simulate_regression_params()
    # booted_slopes = m.simulate_regression_params(boot=True)
    # plot_data = pd.concat(
    #     (
    #         pd.DataFrame(normal_approx, columns=["slope estimate"]).assign(
    #             style="mle asymptotics"
    #         ),
    #         pd.DataFrame(booted_slopes, columns=["slope estimate"]).assign(
    #             style="boostrapped"
    #         ),
    #     )
    # )
    # f, ax = plt.subplots(1, 1, figsize=(6, 4))
    # sns.kdeplot(plot_data, ax=ax, hue="style", x="slope estimate")
    # plt.draw()

    # conditional = m.simulate_predictions()
    # unconditional = m.simulate_predictions(conditional=False)
    # unconditional_booted = m.simulate_predictions(boot=True, conditional=False)

    # plot_data = pd.concat(
    #     (
    #         pd.DataFrame(conditional, columns=["conditional"]),
    #         pd.DataFrame(unconditional, columns=["unconditional"]),
    #         pd.DataFrame(unconditional_booted, columns=["unconditional_booted"]),
    #     ),
    #     axis=1,
    # ).melt(var_name="simulation_type", value_name="replicate")
    # ax = sns.kdeplot(data=plot_data, hue="simulation_type", x="replicate")
    # ax.set_xscale("log")
    # plt.show()


# d = pd.read_csv("/home/phyto/CoccoData/data/unprocessed/poulton2024.csv")
# Y_train = d['pg poc']
# X_train = d['volume']
# X_predict = [100, 99, 98, 101, 102, 100]*1000
# m = regression_simulation(X_train, X_predict, Y_train)
# test1 = m.simulate_predictions()
# sns.histplot(x=test1)
# plt.show()

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


# print("fin")

import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
import glob, os, warnings
from joblib import Parallel, delayed
from scipy import stats
plt.rcParams.update({"font.size": 14})
from functions import rename_synonyms, bayes_bootstrap
from yaml import load, Loader
from collections import namedtuple
import yaml


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

    def __init__(self, X_train, Y_train, 
                 link= sm.families.links.Log()):

        self.link = link
        self.X_train = X_train
        self.Y_train = Y_train
        
        allometric_model = sm.GLM(
            self.Y_train,
            self.X_train,
            family=sm.families.Gamma(link=self.link)       
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
                        family=sm.families.Gamma(link=self.link),
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
        if isinstance(X_predict, pd.DataFrame):
            ncol_X_pred = len(X_predict.columns)
        else:
            ncol_X_pred = 1
        X_predict = np.atleast_2d(X_predict).reshape(-1, ncol_X_pred)

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

    def plot_fit(self, x, y, boot=False, n_replicates=1, 
                    figsize=(12, 8),
                    log_trans = True,
                    ylab = "Carbon content (pg C, log10)",
                    title = ""):

        yhat = pd.DataFrame(
            {
                "yhat": self.simulate_predictions(x,
                    boot=boot, n_replicates=n_replicates
                ).flatten()
            }
        )

        d = pd.DataFrame({'x': x, 
                          'y': y,
                          'yhat': yhat['yhat']})
        
        fig, axs = plt.subplots(1, 2, figsize=figsize)



        #plotting axs[0]:
        #points:
        sns.scatterplot(x=d["x"], y=d["y"],  ax=axs[0])
        #slope:
        slope = self.results.params['volume']

        # Define the line function
        def line(x, slope, intercept):
            return slope * x + intercept
        # # Generate a range of x values for the line
        x_range = np.linspace(min(d["x"]), max(d["x"]), 100)
        y_range = line(x_range, slope, 0)
        axs[0].plot(x_range, y_range, color='steelblue')


        #plotting axs[1]
        sns.scatterplot(x=d["yhat"], y=d["y"], ax=axs[1])
        #fit log model for plotting:
        model = sm.OLS(np.log(d['y']), np.log(d['yhat'])).fit()
        slope = model.params[0]
        axs[1].axline(slope=slope, xy1=[np.min(y), np.min(y)], color="blue")
        #1:1 line
        axs[1].axline(slope=1, xy1=[np.min(y), np.min(y)], color="black", linestyle="dashed")


        if log_trans == True:
            axs[0].set_yscale("log")
            axs[0].set_xscale("log")
            axs[1].set_yscale("log")
            axs[1].set_xscale("log")


        axs[0].set_title("Allometric scaling")
        axs[0].set_ylabel(ylab)
        axs[0].set_xlabel("Cell size (um3, log10)")

        axs[1].set_title("Observed vs Predicted")
        axs[1].set_ylabel("Observed values (pg C, log10)")
        axs[1].set_xlabel("Predicted values (pg C, log10)")

        fig.suptitle(title, size=24, weight="bold")

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

        yhat = pd.DataFrame(
            {
                "yhat": self.simulate_predictions(x,
                    boot=boot, n_replicates=n_replicates
                ).flatten()
            }
        )

        d = pd.concat([x, y, yhat], axis=1)
        d.dropna(inplace=True)
        print(d)

        fig, axs = plt.subplots(2, 2, figsize=figsize)

        if log_trans == True:
            d['volume'] = np.log10(d['volume'])
            d['pg pic'] = np.log10(d['pg pic'])
            d['yhat'] = np.log10(d['yhat'])


        sns.regplot(
            data=d[d["phase_HET"] == 1], x="volume", y="pg pic", ax=axs[0, 0], 
            ci=None, line_kws={"color": "navy"}
        )
        sns.regplot(
            data=d[d["phase_HET"] == 1], y="pg pic", x="yhat", ax=axs[0, 1], 
            ci=None, line_kws={"color": "navy"}
        )

        sns.regplot(
            data=d[d["phase_HOL"] == 1], x="volume", y="pg pic", ax=axs[1, 0], 
            ci=None, color="firebrick"
        )

        sns.regplot(
            data=d[d["phase_HOL"] == 1], y="pg pic", x="yhat", ax=axs[1, 1], 
            ci=None, color="firebrick"
        )


        axs[0, 1].plot(
            d[d["phase_HET"] == 1]["pg pic"], d[d["phase_HET"] == 1]["pg pic"],
            color="black"#, linestyle="dotted"
        )
        axs[1, 1].plot(
            d[d["phase_HET"] == 1]["pg pic"], d[d["phase_HET"] == 1]["pg pic"],
            color="black"#, linestyle="dotted"
        )

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

        axs[0, 1].set_xlim(0, 
                           np.max(d[d["phase_HET"] == 1]['pg pic']))
        
        axs[0, 1].set_ylim(0, 
                           np.max(d[d["phase_HET"] == 1]['yhat']))

        axs[1, 1].set_xlim(0, 
                           np.max(d[d["phase_HET"] == 1]['pg pic']))
        
        axs[1, 1].set_ylim(0, 
                           np.max(d[d["phase_HET"] == 1]['yhat']))
               
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

        fig.suptitle(title)
        plt.tight_layout()
        plt.show()


class library:

    def __init__(self, phases, family, sizes, pic, poc, species_list):
        self.phases_path = phases
        self.family_path = family
        self.sizes_path = sizes
        self.pic_path = pic
        self.poc_path = poc
        self.Study = namedtuple("Study", ['id', 'mean', 'sd', 'method', 'species', 'genera', 'family', 'phase', 'alternate_phase', 'measurement'])

        def def_grouping():

            with open(self.phases_path, 'r') as f:
                phases = load(f, Loader=Loader)

            with open(self.family_path, 'r') as f:
                families = load(f, Loader=Loader)

            d = pd.DataFrame.from_dict(phases, orient='index')
            d = d.rename_axis("species").reset_index()
            d['genera'] = d['species'].str.split(" ").str[0]

            inverse = {}
            for k,v in families.items():
                for x in v:
                    inverse.setdefault(x, []).append(k)

            df = pd.DataFrame.from_dict(inverse, orient='index')
            df = df.rename_axis("genera").reset_index()
            df = df.rename(columns={0: "family"})
            d = pd.merge(d, df, on='genera', how="outer")
            d = d.where(pd.notnull(d), None)

            library = list(d.itertuples(name='species', index=False))

            return(d)

        groups = def_grouping()
        self.species_list = species_list

        def import_data(path):

            all_files = glob.glob(os.path.join(path, "*.csv"))

            d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
            d = d.fillna(0)

            #rename synonyms and typos:
            d = rename_synonyms(d, remove_duplicate=False, check_synonyms = True)

            return(d)

        d = import_data(self.sizes_path)

        d_pic = import_data(self.pic_path)
        d_poc = import_data(self.poc_path)

        list_of_studies = np.concatenate([d['reference'].unique(), 
                                          d_pic['reference'].unique(),
                                          d_poc['reference'].unique()])

        self.list_of_studies = np.unique(list_of_studies)

        def size_and_method(d, study, species):

            d = d[(d['species'] == species) & (d['reference']==study)]

            keys_to_extract = ['mean', 'sd', 'method']
            extracted_data = {}

            for key in keys_to_extract:
                try:
                    extracted_data[key] = d.get(key, None).item()
                except:
                    extracted_data[key] = None
            return extracted_data['mean'], extracted_data['sd'], extracted_data['method']


        def extract_c(d, study, species):
            d = d[(d['species'] == species) & (d['reference']==study)]
            keys_to_extract = ['mean', 'sd', 'method']
            extracted_data = {}
            for key in keys_to_extract:
                try:
                    extracted_data[key] = d.get(key, None).item()
                except:
                    extracted_data[key] = None
            return extracted_data['mean'], extracted_data['sd'], extracted_data['method']

        def classification(groups, species):

            groups = groups[groups['species'] == species]
            keys_to_extract = ['genera', 'family', 'phase', 'alternate_phase']
            extracted_data = {}

            for key in keys_to_extract:
                try:
                    extracted_data[key] = groups.get(key, None).item()
                except:
                    extracted_data[key] = None

            return extracted_data['genera'], extracted_data['family'], extracted_data['phase'], extracted_data['alternate_phase']

        def fill_namedtuple(groups, species_list):
            studies = []
            for id in self.list_of_studies:
                for species in species_list:
                    genera, family, phase, alternate_phase = classification(groups, species)
                    studies.append(self.Study(id, *size_and_method(d, id, species), species, genera, family, 
                                              phase, alternate_phase, measurement="volume"))
                    studies.append(self.Study(id, *extract_c(d_pic, id, species), species, genera, family, 
                                              phase, alternate_phase, measurement="pic"))
                    studies.append(self.Study(id, *extract_c(d_poc, id, species), species, genera, family, 
                                              phase, alternate_phase, measurement="poc"))
            return(studies)

        self.library = fill_namedtuple(groups, species_list)
    
    def return_ntpl(self):
        return(self.library)

    def return_species_list(self):
        return(self.species_list)
    
    def return_HOL_list(self):
        species_library =  [t for t in self.library  if t.phase == "HOL"]
        return([t.species for t in species_library])
    
    def export_yml(self, path):
        spp_list = []
        #sort species list alphabetically:
        for i in range(len(self.species_list)):

            name = self.species_list[i]
            species_library =  [t for t in self.library  if t.species == name]
            species_library_size =  [t for t in self.library  if (t.species == name and
                                                                  t.measurement == "volume")]
            species_library_pic =  [t for t in self.library  if (t.species == name and
                                                                  t.measurement == "pic")]
            species_library_poc =  [t for t in self.library  if (t.species == name and
                                                                  t.measurement == "poc")]
            
            sizes = {study.id:study._asdict() for study in species_library_size }
            for id in self.list_of_studies:
                try:
                    del sizes[id]['id']
                    del sizes[id]['species']
                    del sizes[id]['genera']
                    del sizes[id]['family']
                    del sizes[id]['phase']
                    del sizes[id]['alternate_phase']
                    del sizes[id]['measurement']

                except:
                    None


            #d = asdict(dc, dict_factory=lambda x: {k: v for (k, v) in x if v is not None})

            pic_values = {study.id:study._asdict() for study in species_library_pic }
            for id in self.list_of_studies:
                try:
                    del pic_values[id]['id']
                    del pic_values[id]['species']
                    del pic_values[id]['genera']
                    del pic_values[id]['family']
                    del pic_values[id]['phase']
                    del pic_values[id]['alternate_phase']
                    del pic_values[id]['measurement']

                except:
                    None


            poc_values = {study.id:study._asdict() for study in species_library_poc }
            for id in self.list_of_studies:
                try:
                    del poc_values[id]['id']
                    del poc_values[id]['species']
                    del poc_values[id]['genera']
                    del poc_values[id]['family']
                    del poc_values[id]['phase']
                    del poc_values[id]['alternate_phase']
                    del poc_values[id]['measurement']

                except:
                    None

            species =  {name: {
                    'genera': species_library[0].genera,
                    'family': species_library[0].family,
                'phase': species_library[0].phase,
                    'alternate_phase': species_library[0].alternate_phase,
                    'size' : sizes,
                    'pic' : pic_values
                }}
            spp_list.append(species)

        with open(path, 'w') as outfile:
            yaml.dump(spp_list, outfile, default_flow_style=False)

        print("exported yml to: " + str(path))



class pipeline:

    def __init__(self, root= '../CASCADE/data/', n = 1000):

        d = pd.read_csv(root + "output/abundant_species.csv")
        self.species_list = d['species']
        self.root = root

        m = library(root + 'classification/phases.yml',
                    root + 'classification/family.yml',
                    root + "sizes/",
                    root + "pic/",
                    root + "poc/",           
                    self.species_list)
        
        m.export_yml(root + "/output/library.yml")

        #create a list of HOL species for future use:
        self.HOL_list = m.return_HOL_list()

        #return the named tuple so we can use it:
        self.ntpl = m.return_ntpl()
        self.ntpl_size =  [t for t in self.ntpl  if t.measurement == "volume"]
        self.ntpl_pic =  [t for t in self.ntpl  if t.measurement == "pic"]
        self.ntpl_poc =  [t for t in self.ntpl  if t.measurement == "poc"]

        self.n = n 
    def return_species_list(self):
        return(self.species_list)

    def check_HOL(self, name):
        """
        returns name of HET phase 
        if size data of HOL is undefined
        and if alternative phase is defined
        """
        species_library =  [t for t in self.ntpl  if t.species == name]
        if all(d.mean is None for d in species_library): 
            if (any(d.phase == "HOL" for d in species_library)):
                if (any(d.alternate_phase is not None for d in species_library)):
                    return(species_library[0].alternate_phase)
        else:
            return(None)

    def resample_cv(self, size_library):
        print("estimating cross study cv for spp")
        #estimate mean SD across studies by bootstrapping
        sd_library = np.array([ x.sd for x in size_library if (x.sd is not None) and (x.mean is not None)])
        mean_library = np.array([x.mean for x in size_library if (x.sd is not None) and (x.mean is not None) ])
        #drop entries where sd values are None
        sd_library = sd_library[sd_library != np.array(None)]
        mean_library = mean_library[sd_library > 0]
        sd_library = sd_library[sd_library > 0]
        #convert sd to cv by dividing by mean:
        cv_library = sd_library/mean_library

        #bootstrap cv estimates
        if cv_library.any():
            cv_estimate = bayes_bootstrap(cv_library, self.n)
        else:
            cv_estimate = [0]*self.n

        return(cv_estimate)
        
    def resample_size(self, spp_name):

        print("resampling size for spp: " + spp_name)

        if self.check_HOL(spp_name) !=None:
            HET_name = self.check_HOL(spp_name)
            size_library =  [t for t in self.ntpl  if t.species == HET_name]
            print(HET_name)
        else:
            size_library =  [t for t in self.ntpl  if t.species == spp_name]

        #simulate SDs for studies which do not have a SD:
        all_studies_cv = self.resample_cv(size_library)
        all_species_cv = self.resample_cv([t for t in self.ntpl])

        #simulate size distributions for each study which has data: 
        estimate = []

        for i in range(len(size_library)):
            #for every study:
            if (size_library[i].mean is not None) and ((size_library[i].sd is None) or (size_library[i].sd ==0)):
                size_estimate = []
                if np.mean(all_studies_cv)  > 0:
                    for j in range(self.n):
                        random_sd =  all_studies_cv[j]*size_library[i].mean
                        size_estimate_n = np.random.normal(size_library[i].mean, random_sd, 1)
                        size_estimate.extend(size_estimate_n)
                else:
                    print("group sd unknown, using cross species sd instead (size)")
                    for j in range(self.n):
                        #if 0 or NA:
                        random_sd =  all_species_cv[j]*size_library[i].mean
                        size_estimate_n = np.random.normal(size_library[i].mean, random_sd, 1)
                        size_estimate.extend(size_estimate_n)

            elif (size_library[i].mean is None) and (size_library[i].sd is None):
                size_estimate = []
            else:
                size_estimate = np.random.normal(size_library[i].mean, size_library[i].sd, self.n)
            #    print(size_library[i].id)
            size_estimate = [x for x in size_estimate if x > 0 ]
            estimate.extend(size_estimate)


        #bayes_bootstrap to merge estimates
        estimate = bayes_bootstrap(estimate)
        return(estimate)


    def c_regression(self, name, measurement, X_predict):
        #function to estimate carbon based on size distribution
        if measurement == "pic":
            d = pd.read_csv(self.root + "/allometry/sheward2024.csv")

            phase = [t.phase for t in self.ntpl  if t.species == name][0]

            d = pd.get_dummies(d, columns=["phase"], dtype=float)
            X_train = d[["volume", "phase_HET", "phase_HOL"]].astype(float)
            Y_train = d['pg pic']

            if phase=='HET' or phase=='CER' or phase=='NANO':
                X_predict_oh = pd.DataFrame(
                    {
                        "volume": X_predict,
                        "phase_HET": [int(1)] * len(X_predict),
                        "phase_HOL": [int(0)] * len(X_predict),
                    }
                )
            elif phase=='HOL' or phase=='POL':
                X_predict_oh = pd.DataFrame(
                    {
                        "volume": X_predict,
                        "phase_HET": [int(0)] * len(X_predict),
                        "phase_HOL": [int(1)] * len(X_predict),
                    }
                )
            else:
                raise ValueError("phase not defined: " + phase)         

            r = regression_simulation(X_train, Y_train=Y_train, link= sm.families.links.Identity())
            c_simulation = r.simulate_predictions(X_predict_oh, boot=True, n_replicates=1)

        elif measurement == "poc":
            d = pd.read_csv(self.root + "/allometry/poulton2024.csv")
            Y_train = d['pg poc']
            X_train = d['volume']

            #X_train = np.log(X_train)
            #X_predict = np.log(X_predict)

            #r = regression_simulation(X_train, Y_train)
            r = regression_simulation(X_train, Y_train=Y_train, link= sm.families.links.Identity())
            c_simulation = r.simulate_predictions(X_predict, boot=True, n_replicates=1)

        else:
            raise ValueError("measurement undefined, should be 'pic' or 'poc'")


        c_simulation = c_simulation.flatten()

        return(c_simulation)


        # try:
        #     c_simulation = c_simulation.str[0]
        #     print(c_simulation)
        # except:
        #     None

        return(c_simulation)


        
    def resample_carbon_species(self, spp_name, size_simulation, measurement):

        library =  [t for t in self.ntpl  if (t.species == spp_name and
                                            t.measurement == measurement)]

        print("estimating " + measurement + " for spp " + spp_name)

        #simulate size distributions for each study which has data: 
        estimate = []

        #drop empty studies from library
        library =  [t for t in library  if t.mean != None]

        if len(library)>0:
            try:
                #simulate SDs of all species to be used for studies which do not have a SD:
                all_studies_cv = self.resample_cv(library)
                all_species_cv = self.resample_cv([t for t in self.ntpl])

            #if direct measurements exist use them:
                print("estimating "+ measurement + " based on resampling")
                for i in range(len(library)):
                    #if mean is known:
                    if (library[i].sd is None) or (library[i].sd ==0):
                        print("mean known for:" + spp_name)
                        #if sd is not known
                        if np.mean(all_studies_cv)  > 0:
                            #if species sd is known
                            for j in range(self.n):
                                random_sd =  all_studies_cv[j]*library[i].mean
                                estimate.extend(np.random.normal(library[i].mean, random_sd, 1))
                        else:
                            #if species sd is unknown
                            print("group sd unknown, using cross species sd instead")
                            for j in range(self.n):
                                #if 0 or NA:
                                random_sd =  all_species_cv[j]*library[i].mean
                                estimate.extend(np.random.normal(library[i].mean, random_sd, 1))
                    else:
                        #if sd is known, use study sd
                        estimate.extend( np.random.normal(library[i].mean, library[i].sd, self.n))
            except:
                #if resample of sd fails because there are no SDs
                #estimate carbon based on size:
                print("estimating "+ measurement + " based on GLM")
                estimate.extend(self.c_regression(spp_name, measurement, size_simulation))   

            estimate = bayes_bootstrap(estimate)
                       
        else:
        #otherwise estimate carbon based on size:
            print("estimating "+ measurement + " based on GLM")
            estimate.extend(self.c_regression(spp_name, measurement, size_simulation))    

        
        return(estimate)



    def resample_measurements(self, species_name):

        size_simulation = self.resample_size(species_name) 
        diameter_simulation = (6*np.asarray(size_simulation)/np.pi)**(1/3)

        pic_simulation = self.resample_carbon_species(species_name, size_simulation, "pic")
        poc_simulation = self.resample_carbon_species(species_name, size_simulation, "poc")

        d_vol = pd.DataFrame({'species': species_name, 'value': size_simulation, 'variable':'volume (um3)'})
        d_poc = pd.DataFrame({'species': species_name, 'value': poc_simulation, 'variable':'pg poc'})
        d_dia = pd.DataFrame({'species': species_name, 'value': diameter_simulation, 'variable':'diameter (um)'})
        d_pic = pd.DataFrame({'species': species_name, 'value': pic_simulation, 'variable':'pg pic'})

        d = pd.concat([d_dia, d_vol, d_poc, d_pic])
        #d = pd.concat([d_dia, d_vol, d_poc])

        return(d)
    
    def export_all(self, export_path, species = None, loop_range = None):

        species_list = self.return_species_list()

        if loop_range==None:
            loop_range = range(0, len(species_list))

        if species == None:
            estimates = []
            for i in loop_range: #
                print(species_list[i])
                estimates.append(self.resample_measurements(species_list[i]))
                print("finished estimating species #" + str(i+1) + " out of " + str(len(species_list)) + " species")
            estimates = pd.concat(estimates)
            estimates.to_csv(export_path, index=False)

        else:
            t = self.resample_measurements(species)
            t.to_csv(export_path, index=False)

class merge_abundances():

    def __init__(self, import_path, export_path, synonym_path):
        self.export_path = export_path
        all_files = glob.glob(os.path.join(import_path, "*.csv"))

        d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)

        d.reset_index(drop=False, inplace=True)

        d['Month'] = d['Month'].astype('int64')
        d['Year'] = d['Year'].astype('int64')

        d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'], inplace=True)

        d = d.replace({0:pd.NA})

        d.rename(columns=lambda x: x.strip(), inplace=True)

        with open(synonym_path, 'r') as f:
            groupings = load(f, Loader=Loader)


        dict = {species:k
            for k, v in groupings.items()
            for species in v}
        try:
            d = d.drop(columns=['Phaeocystis pouchetii'])
        except:
            None
        try:
            d = d.drop(columns=['Thoracosphaera heimii'])
        except:
            None

        d = d[d.sum(axis=1)>0]

        d = (d.rename(columns=dict)
            .groupby(level=0, axis=1, dropna=False)).sum( min_count=1).reset_index()

        self.references = d['Reference'].unique()

        print(self.references)
        #d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).agg('mean')
        #d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).mean()

        #species = d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).columns
        #d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Month', 'Year', 'Reference', 'Method'])[species].agg('mean') #.reset_index()
        d.set_index(['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method'], inplace=True)
        print("groupby")
        #print(d.reset_index()['Reference'].unique())


        d = d.drop(columns=['index'])

        print("drop index")
        #print(d.reset_index()['Reference'].unique())

        #check if final species are in species library
        species_library = {v: k for k, v in dict.items()}

        species_observed = d.columns

        if (set(species_observed).issubset(species_library)):
            print("all species are defined")
        else:
            raise ValueError("undefined species:" + str(set(species_observed).difference(species_library)))

        #drop any columns that contain "undefined"
        d = d[d.columns.drop(list(d.filter(regex='undefined')))]

        try: #make new dir if needed
            os.makedirs(self.export_path)
        except:
            None

        print("dropped columns that contained undefined")
        #print(d.reset_index()['Reference'].unique())

        counts =pd.DataFrame({'count': np.count_nonzero(d.fillna(0), axis=0), 'species': d.columns})
        counts = counts[counts['species']!="Reticulofenestra sessilis"]
        filter = counts['species'].str.contains('undefined')
        counts = counts[~filter]
        counts = counts[counts['count']>0]

        non_zero_spp = counts[counts['count']>0]['species']

        d = d[non_zero_spp]
        #print(d['Reference'].unique())

        self.counts = counts
        self.abundant = counts[counts['count']>=20]
        self.rare = counts[counts['count']<20]

        self.non_zero_spp = non_zero_spp

        self.d = d


    def export_obs_counts(self): 
        self.counts.sort_values(by=["count"], ascending=False).to_csv(self.export_path + "counts.csv", index=False )
        self.abundant.sort_values(by=["count"], ascending=False).to_csv(self.export_path + "abundant_species.csv", index=False )
        self.rare.sort_values(by=["count"], ascending=False).to_csv(self.export_path + "rare_species.csv", index=False )


    def export_ungridded_abundances(self):
        df = self.d #[self.non_zero_spp]
    
        df.to_csv(self.export_path +"ungridded_abundances.csv")
        print("ungridded abundances exported to: " + self.export_path +"ungridded_abundances.csv")

    def gridding(self):
        d = self.d[self.abundant['species']]
        d = d.reset_index()
        depth_bins = np.linspace(-1, 300, 62).astype(np.int64) 
        depth_labels = np.linspace(0, 300, 61).astype(np.int64) 
        d = d[d["Depth"] >= 0]
        d = d[d["Depth"] <= 301]

        d['Depth'] = pd.cut(d['Depth'], bins=depth_bins, labels=depth_labels).astype(np.int64) 

        lat_bins = np.linspace(-90, 90, 181)
        lat_labels = np.linspace(-90, 89, 180)
        d['Latitude'] = pd.cut(d['Latitude'].astype(np.float64), bins=lat_bins, labels=lat_labels).astype(np.float64) 

        lon_bins = np.linspace(-180, 180, 361)
        lon_labels = np.linspace(-180, 179, 360)
        d['Longitude'] = pd.cut(d['Longitude'].astype(np.float64), bins=lon_bins, labels=lon_labels).astype(np.float64) 

        species = d.drop(columns=['Latitude', 'Longitude', 'Depth', 'Day', 'Month', 'Year', 'Reference', 'Method']).columns

        d = d.groupby(['Latitude', 'Longitude', 'Depth', 'Month', 'Year'])[species].agg('mean').reset_index()
        d = d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year'])
    
        d = d.reset_index()

        return(d)


    def export_gridded_abundances(self):
        d = self.gridding()

        d.to_csv(self.export_path +"gridded_abundances.csv", index=False)
        
        print("gridded abundances exported to: " + self.export_path +"gridded_abundances.csv")



    def print_stats(self):

        df = self.d[self.non_zero_spp]
        print("concatenated abundances: " +str(len(df)))

        t = pd.melt(df)
        t = t[t['value']>0]

        print("concatenated observations: " +str(len(t)))

        d = self.gridding()
        print("gridded abundances: " +str(len(d)))

        d.set_index(['Latitude', 'Longitude', 'Depth', 'Month', 'Year'])
        t = pd.melt(d)
        t = t[t['value']>0]

        print("gridded observations: " +str(len(t)))

        
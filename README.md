# CASCADE

CASCADE is a global dataset for 139 extant coccolithophore taxonomic units. CASCADE includes a trait database (size and cellular organic and inorganic carbon contents) and taxonomic-specific global spatiotemporal distributions (Latitude/Longitude/Depth/Month/Year) of coccolithophore abundance and organic and inorganic carbon stocks. CASCADE covers all ocean basins over the upper 275 meters, spans the years 1964-2019 and includes 33,119 taxonomic-specific abundance observations. Within CASCADE, we characterise the underlying uncertainties due to measurement errors by propagating error estimates between the different studies. 

## Installing the dependencies (w/ conda):

Install the dependencies in a new environment: 

``` conda env create -f cascade_save_path/environment.yml ``` 

This will create a conda environment called "cascade-env" which can then be loaded in your Jupyter notebook environment.


## Reproducing the Zenodo CASCADE data archive:

To create the final data set simply run the notebooks provided in ```./notebooks/```.

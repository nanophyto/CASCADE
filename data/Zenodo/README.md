This is the Zenodo data archive for the Coccolithophore Abundance, Size, Carbon And Distribution Estimates (CASCADE) dataset. CASCADE is a global dataset for 139 extant coccolithophore taxonomic units. CASCADE includes a trait database (size and cellular organic and inorganic carbon contents) and taxonomic-specific global spatiotemporal distributions (Latitude/Longitude/Depth/Month/Year) of coccolithophore abundance and organic and inorganic carbon stocks. CASCADE covers all ocean basins over the upper 275 meters, spans the years 1964-2019 and includes 31,875 taxonomic-specific abundance observations. Within CASCADE, we characterise the underlying uncertainties due to measurement errors by propagating error estimates between the different studies. 

Full details of the data set are provided in the associated Scientific Data manuscript. The repository contains five main folders: 1) "Classification", which contains YAML files with synonyms, family-level classifications, and life cycle phase associations and definitions; 2) "Concatenated literature", which contains the merged datasets of size, PIC and POC and which were corrected for taxonomic unit synonyms; 3) "Resampled cellular datasets", which contains the resampled datasets of size, PIC and POC in long format as well as a summary table; 4) "Gridded data sets", which contains gridded datasets of abundance, PIC and POC; 5) "Species list", which contains spreadsheets of the "common" (>20 obs) and "rare" (<20 obs) species and their number of observations. 

The CASCADE data set can be easily reproduced using the scripts and data provided in the associated github repository: https://github.com/nanophyto/CASCADE.

```bash

CASCADE
├── classification
|   └── family.yml
|   └── phases.yml
|   └── synonyms.yml
├── concatenated_literature
|   └── literature_abundances.csv
|   └── literature_PIC.csv
|   └── literature_POC.csv
|   └── literature_sizes.csv
├── resampled_cellular_datasets
|   └── resampled_diameter.csv 
|   └── resampled_PIC.csv 
|   └── resampled_POC.csv 
|   └── resampled_volume.csv 
|   └── summary_table.csv
├── gridded_datasets
|   └── gridded_abundances.scv
|   └── gridded_pic.scv
|   └── gridded_poc.scv
├── species_lists
|   └── common_species.csv 
|   └── rare_species.csv 
└── README.md

```
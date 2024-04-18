# UrbanFoodChain
A repository that contains the data and code for:

Larson, R.N., M. Fidino, and H. A. Sander. Urban deer mouse abundance is more strongly correlated with predator occupancy than land cover. *Landscape and Urban Planning*. (In review)


This `README` file includes information on the various scripts and datasets used for this analysis. Not every data source is saved in this repository (e.g., GIS data). The manuscript contains citations for where to access the geospatial data.


<h2>Scripts</h2> </div>

**./Mouse_Model.R:** R script that cleans and processes the data & runs the population model for the deer mice

**./Predator_Models.R:** R script that cleans and processes the data, runs the multi-species occupancy model, and performs post-model calculations (e.g., average occupancy probability) for the predators (coyote, red fox, domestic cats, mink)

**./graphing.R:** R script that graphs the mouse model results (Figures 2 & 3 in the manuscript, and Supplemental Figures 1 & 2)

<h2>Folders</h2>
<h3>data</h3>

There is 1 subfolder and 4 files in this folder

**./data/PERO.csv:** The counts of deer mice captured for each night of trapping

Data in this file take the following format. Cells that contain an 'NA' instead of a number mean traps were not set at that site for those nights (i.e., the site was not sampled):
| Column       | Type    | Description                                                                                                      |
| ------------ | ------- | ---------------------------------------------------------------------------------------------------------------- |
| Site         | numeric | The identity/name of each trapping plot.                                                                         |
| SP21_1 | numeric | The counts of individuals of that species for the 1st night of trapping during spring 2021                                                  |
| SP21_2 | numeric | The counts of individuals of that species for the 2nd night of trapping during spring 2021                                                  |
| SP21_3 | numeric | The counts of individuals of that species for the 3rd night of trapping during spring 2021                                                  |
| SU21_1 | numeric | The counts of individuals of that species for the 1st night of trapping during summer 2021                                                  |
| ... | numeric | ...                                                |
|y_x  | numeric | The counts of individuals of that species for the x night of trapping during y sample period                                     |

**./data/predatorOccu.csv:** Detection histories of the predaors included in the predator model

Data in this file take the following format:
| Column       | Type    | Description                                                                                                      |
| ------------ | ------- | ---------------------------------------------------------------------------------------------------------------- |
| Species         | string | The species (coyote, cat, fox, or mink)                                                                       |
| Sample_Period | numeric | The sampling period, from 1 (Fall 2019) to 13 (Fall 2022) |
| Site | string | The name of the camera site                                             |
| Y | numeric | The number of detections during the sampling period of the species from the 'Species' column                                                |
| J  | numeric | The number of days the camera was active during the sampling period                                     |

**./data/siteCovs_Pred.csv:** The site covariates for each of the 41 camera trap locations. All covariates were measured from a 500m buffer surrounding the camera's location.

| Column    | Type    | Description                                                                                                        |
| --------- | ------- | ------------------------------------------------------------------------------------------------------------------ |
| Site  | string | The name of the camera site (matches name in the predatorOccu data table).                                                                           |
| Imperv    | numeric | The average percent cover of impervious surfaces in the 500m buffer                                                                     |
| Turfgrass     | numeric | The percent cover of turfgrass (NLCD categories "Developed, Open Space"; and "Developed, Low Intensity" combined) in the 500m buffer |
| Forest  | numeric | The percent cover of forest (NLCD categories "Coniferous Forest", "Deciduous Forest", "Mixed Forest", and "Woody Wetlands" combined) in the 500m buffer |
| Prairie  | numeric | The percent cover of grasslands (NLCD categories "Herbacous Vegetation", "Hay/Pasture", and "Emergent Herbacous Wetlands" combined) in the 500m buffer                                                                      |
| Crop   | numeric | The percent cover of crops (NLCD category "Cultivated Crops") in the 500m buffer |
| ResUnits | numeric | The number of residential housing units in the 500m buffer |
| Dist_to_Wat    | numeric | The distance from the camera location to the nearest water body (stream, lake, or river) in meters |

**./data/siteCovs_Rod.csv:** The site covariates for each of the 45 small mammal trapping plots.

| Column    | Type    | Description                                                                                                        |
| --------- | ------- | ------------------------------------------------------------------------------------------------------------------ |
| site  | numeric | The identity of each trapping plot.                                                                           |
| canopy    | numeric | The average tree canopy closure on each plot                                                                       |
| shrub     | numeric | The average percent cover of vegetation between 76 - 500 cm in height on each plot |
| tallHerb  | numeric | The average percent cover of vegetation between 0 - 75 cm in height on each plot |
| humanMod  | numeric | Sum of the 'imperv' and 'turfgrass' columns                                                                        |
| imperv    | numeric | The average percent cover of impervious surfaces on each plot |
| turfgrass | numeric | The average percent cover of turfgrass on each plot |
| contag    | numeric | The contagion index of each plot (i.e., a measure of habitat fragmentation). |

<h3>data/obsVars</h3>

There are 3 files in this subfolder

**./data/obsVars/jDate.csv:** The Julian date of each night of small mammal trapping

**./data/obsVars/Moon.csv:** The moon illumination (proportion full) on each night of trapping

**./data/obsVars/Effort.csv:** The number of available traps (i.e., traps that remained undisturbed) on each night of trapping

All files follow this format, with the jDate file as an example. Just sub out 'jDate' for 'moon' with the moon illumination data and 'effort' for the trap effort data:
| Column       | Type    | Description                                                                                                      |
| ------------ | ------- | ---------------------------------------------------------------------------------------------------------------- |
| Site         | numeric | The identity/name of each trapping plot.                                                                         |
| jDate_SP21_1 | numeric | The Julian Date of the 1st night of trapping during spring 2021                                                  |
| jDate_SP21_2 | numeric | The Julian Date of the 2nd night of trapping during spring 2021                                                  |
| jDate_SP21_3 | numeric | The Julian Date of the 3rd night of trapping during spring 2021                                                  |
| jDate_SU21_1 | numeric | The Julian Date of the 1st night of trapping during summer 2021                                                  |
| ... | numeric | ...                                                |
| jDate_y_x  | numeric | The Julian Date of the x night of trapping during y sampling period                                                 |

<h3>JAGS</h3>

There is 1 file in this folder, `./JAGS/dynamicCommunityModel.R`, which is the Bayesian community abundance model code to pass to JAGS.

<h3>functions</h3>

There are 4 files in this folder, all utility functions that automate or declutter some of the R code

**.functions/logit2prob.R:** Function script to convert a logit value to a probability

**.functions/split_mcmc.R:** Function to split a model's MCMC matrix into a list of named objects, one for every parameter. Makes graphing results much easier. Credit for this code goes to [@mfidino](https://github.com/mfidino) [(see his blog post here)](https://masonfidino.com/split_mcmc/)

**.functions/wide_to_stacked.R:** Function to convert a wide-format data frame (i.e., one site per row with one column for each observation) into a stacked format data frame (one observation per row, with sites/seasons/etc. "stacked" on top of each other). Modified from code in [a vignette for the `umbs` R package](github.com/kenkellner/umbs/blob/master/vignettes/random-effects.Rmd).

**.functions/wideObs_to_stacked.R:** Function script to convert observation-level covariate data from a wide format to a 'stacked' format. Different from `wide_to_stacked.R` in that it does not add a 'Species' column to the resulting dataframe. Modified from code in [a vignette for the `umbs` R package](github.com/kenkellner/umbs/blob/master/vignettes/random-effects.Rmd).

<h3>landscapes</h3>

There is 1 subfolder and 1 file in this folder. The contents of this folder are for creating the contagion index values for each site (see the manuscript for more details).

**./landscapes/connectivity.R:** R script for calculating habitat fragmentation (i.e., contagion index) for each trapping plot

<h3>landscapes/tifs</h3>

This subfolder contains `.tif` files of the land cover of each small mammal trapping plot. There are 45 files in here, so I'm not going to list them all, but they are named after their site number: e.g., `1.tif`, `124.tif`, etc.

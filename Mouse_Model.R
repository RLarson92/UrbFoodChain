#########################################################
#                                                       #
#           Dynamic Mouse Abundance Model               #
#                                                       #
#           Last Updated: 7 June 2024                   #
#                                                       #
#########################################################
# this model is a population abundance model with fully
# modeled 'recruitment' (births + immigrants) and 
# 'survival' (1-[deaths + emigrants])

#### Data Cleaning & Preparation ####
library(tidyr)
library(dplyr)
# load functions used to clean data
functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}
# read in the mouse count dataset
df <- read.csv("./data/PERO.csv")
dets <- df[ ,2:19] # 2:19 for full dataset, 2:7 for 2-season test of code
# stack the data, long-style
stack <- wide_to_stacked(dets, 3)
# stack the replicates (nights) on top of each other and order the datasheet. 
# Data should be in order by season then site
n_long <- gather(stack, key = "night", value = "N", night1:night3, factor_key = TRUE)
n_long <- n_long[order(
  n_long[,"Season"],
  n_long[,"Site"]), ]
n_long$night <- as.numeric(n_long$night)
n_long$Site <- as.numeric(as.factor(n_long$Site))

# create a blank array that has the correct dimensions
my_array <- array(
  NA,
  dim = c(
    max(n_long$Site),
    max(n_long$night),
    max(n_long$Season)
  )
)
# fill in the array with the information from the long-form data.frame
for(i in 1:nrow(n_long)){
  my_array[
    n_long$Site[i],
    n_long$night[i],
    n_long$Season[i]
  ] <- n_long$N[i]
}
rm(dets, stack)

# Site covariates are a 2-D array [sites, covariates]
site_covs <- read.csv("./data/siteCovs_Rod.csv", stringsAsFactors = FALSE)
siteCovs <- scale(site_covs[ ,3:6])
# these variables are co-linear, so we'll do a PCA to collapse them
pca <- princomp(siteCovs)
pca$loadings
# PC1: more + = canopy/herb/shrub; more - = humanMod
# pull in contagion indices & join. assumes you have run the 'Connectivity.R'
# code in the 'landscapes' folder first
contag <- read.csv("./data/contag.csv", stringsAsFactors = FALSE)
# pull in modeled predator occupancy (this code assumes you have run
# the 'Predator_Models.R" code first)
occuPred <- read.csv(./data/modeledOccu_Predators.csv", stringsAsFactors = FALSE)
Pred <- occuPred[,2:4]
Pred <- cbind(Pred, occuPred[,7:8])
# let's make our covariate data.frame
covs_for_model <- cbind.data.frame("PC1" = pca$scores[,1],
                                   "site" = site_covs[,1],
                                   "camera" = site_covs[,2])
covs_for_model <- left_join(contag, covs_for_model, by="site")
covs_for_model$contag <- scale(covs_for_model$contag)
covs_for_model <- left_join(Pred, covs_for_model, by="camera")
covs_for_model <- drop_na(covs_for_model)
covs_for_model <- covs_for_model[,-1]
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
covs_for_model <- covs_for_model[order(
  covs_for_model[,"site"]
),]
covMatrix <- as.matrix(covs_for_model)
rm(site_covs)

# Observation covariates are 3-D array [sites, nights, seasons]
# We have 3 separate arrays, one for each covariate
moon <- read.csv("./data/obsVars/Moon.csv", stringsAsFactors = FALSE)
moon1 <- as.data.frame(moon[,2:19]) # to 19 for full dataset, to 10 for test
moonStack <- wideObs_to_stacked(moon1, 3) 
moon_long <- gather(moonStack, key = "night", value = "moon", night1:night3, factor_key = TRUE)
moon_long$moon <- as.numeric(scale(moon_long$moon))
moon_long <- moon_long[order(
  moon_long[,"Season"],
  moon_long[,"Site"]),]
moon_long$night <- as.numeric(moon_long$night)
moon_long$Site <- as.numeric(as.factor(moon_long$Site))
moonArray <- array(
  NA,
  dim = c(max(moon_long$Site),
          max(moon_long$night),
          max(moon_long$Season)
  )
)
for(i in 1:nrow(moon_long)){
  moonArray[
    moon_long$Site[i],
    moon_long$night[i],
    moon_long$Season[i]
  ] <- moon_long$moon[i]
}
rm(moon1, moonStack, moon_long)

Date <- read.csv("./UrbFoodChain/data/obsVars/Date.csv", stringsAsFactors = FALSE)
Date1 <- as.data.frame(Date[,2:19]) # to 19 for full dataset, to 10 for test
DateStack <- wideObs_to_stacked(Date1, 3) 
Date_long <- gather(DateStack, key = "night", value = "Date", night1:night3, factor_key = TRUE)
Date_long$Date <- as.numeric(scale(Date_long$Date))
Date_long <- Date_long[order(
  Date_long[,"Season"],
  Date_long[,"Site"]),]
Date_long$night <- as.numeric(Date_long$night)
Date_long$Site <- as.numeric(as.factor(Date_long$Site))
DateArray <- array(
  NA,
  dim = c(max(Date_long$Site),
          max(Date_long$night),
          max(Date_long$Season)
  )
)
for(i in 1:nrow(Date_long)){
  DateArray[
    Date_long$Site[i],
    Date_long$night[i],
    Date_long$Season[i]
  ] <- Date_long$Date[i]
}
rm(Date1, DateStack, Date_long, Date)

effort <- read.csv("./data/obsVars/Effort.csv", stringsAsFactors = FALSE)
effort1 <- as.data.frame(effort[,2:19]) # to 19 for full dataset, to 10 for test
effortStack <- wideObs_to_stacked(effort1, 3) 
effort_long <- gather(effortStack, key = "night", value = "effort", night1:night3, factor_key = TRUE)
effort_long$effort <- as.numeric(scale(effort_long$effort))
effort_long <- effort_long[order(
  effort_long[,"Season"],
  effort_long[,"Site"]),]
effort_long$night <- as.numeric(effort_long$night)
effort_long$Site <- as.numeric(as.factor(effort_long$Site))
effortArray <- array(
  NA,
  dim = c(max(effort_long$Site),
          max(effort_long$night),
          max(effort_long$Season)
  )
)
for(i in 1:nrow(effort_long)){
  effortArray[
    effort_long$Site[i],
    effort_long$night[i],
    effort_long$Season[i]
  ] <- effort_long$effort[i]
}
rm(effort1, effortStack, effort_long, effort, moon)

# finally, coding up a season dummy variable to account for increasing
# mouse populaions as young are born throughout the year
seasonData <- as.factor(c("spring","summer"
                          ,"fall","spring","summer","fall"
                          ))

#### RUN MODEL ####
# Data list for model
data_list <- list(
  nsite = max(n_long$Site),
  PC1 = covMatrix[,7],
  contag = covMatrix[,5],
  cat.mu = covMatrix[,1],
  cat.sd = covMatrix[,2],
  fox.mu = covMatrix[,3],
  fox.sd = covMatrix[,4],
  nNight = max(n_long$night),
  moon = moonArray,
  jDate = jDateArray,
  effort = effortArray,
  y = my_array, 
  nSP = max(n_long$Season), #nSP stands for 'number of sample periods' to distinguish it from 'season' representing meteorological season
  season = seasonData,
  nseason = 3
)
save(data_list, file = "./results/data_list.Rdata") #not really a 'result' but saving there for later use

# initial values for N should be maximum possible counts
N_init = array(10, dim = c(data_list$nsite, data_list$nSP))
# set to NA for all seasons t>1 because we cannot know this (needs to be specified by R + S)
N_init[,2:6] <- NA
# initial values for recruits can also be set to maximum possible counts
R_init = array(10, dim = c(data_list$nsite, data_list$nSP))
# providing null values for t = 1 (i.e., no recruits yet)
R_init[,1] <- NA

# Fit the model
library(runjags)
my_mod <- runjags::run.jags(
  model = "./JAGS/mouseModel.R",
  monitor = c("beta0","beta1",
              "phi0", "phi1", "phi2", "phi3", "phi4", "phi5", "phi6",
              "gamma0", "gamma1", "gamma2",
              "alpha0", "alpha1", "alpha2", "alpha3",
              "cat","fox"
              ),
  data = data_list,
  n.chains = 3,
  inits = rod_inits,
  burnin = 100000,
  sample = 25000,
  adapt = 1000,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
# save the output for later & print the summary to check for convergence, etc.
saveRDS(my_mod, "./results/rod_model.RDS")
summary(my_mod, vars = c("beta0","beta1",
                         "phi0", "phi1", "phi2", "phi3", "phi4", "phi5", "phi6",
                         "gamma0", "gamma1", "gamma2",
                         "alpha0", "alpha1", "alpha2", "alpha3"
))
plot(my_mod)


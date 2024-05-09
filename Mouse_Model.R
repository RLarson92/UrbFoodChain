library(tidyr)

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
siteCovs <- scale(site_covs[ ,2:5])
# these variables are co-linear, so we'll do a PCA to collapse them
pca <- princomp(siteCovs)
pca$loadings
# PC1: more + = canopy/herb/shrub; more - = humanMod
# pull in modeled predator occupancy (this code assumes you have run
# the 'Predator_Models.R" code first
occuPred <- read.csv(./data/modeledOccu_Predators.csv", stringsAsFactors = FALSE)
occuPred$cat.sd <- abs(site_covs$cat.sd)
occuPred$fox.sd <- abs(site_covs$fox.sd)
# let's make our covariate data.frame
covs_for_model <- cbind.data.frame("PC1" = pca$scores[ ,1],
                                   "contag" = scale(site_covs$contag),
                                   "site" = site_covs[ ,1])
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
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

jDate <- read.csv("./data/obsVars/jDate.csv", stringsAsFactors = FALSE)
jDate1 <- as.data.frame(jDate[,2:19]) # to 19 for full dataset, to 10 for test
jDateStack <- wideObs_to_stacked(jDate1, 3) 
jDate_long <- gather(jDateStack, key = "night", value = "jDate", night1:night3, factor_key = TRUE)
jDate_long$jDate <- as.numeric(scale(jDate_long$jDate))
jDate_long <- jDate_long[order(
  jDate_long[,"Season"],
  jDate_long[,"Site"]),]
jDate_long$night <- as.numeric(jDate_long$night)
jDate_long$Site <- as.numeric(as.factor(jDate_long$Site))
jDateArray <- array(
  NA,
  dim = c(max(jDate_long$Site),
          max(jDate_long$night),
          max(jDate_long$Season)
  )
)
for(i in 1:nrow(jDate_long)){
  jDateArray[
    jDate_long$Site[i],
    jDate_long$night[i],
    jDate_long$Season[i]
  ] <- jDate_long$jDate[i]
}
rm(jDate1, jDateStack, jDate_long, jDate)

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

seasonData <- as.factor(c("spring","summer"
                          ,"fall","spring","summer","fall"
                          ))

#### RUN MODEL ####
# Data list for model
data_list <- list(
  nsite = max(n_long$Site),
  PC1 = covMatrix[,2],
  contag = covMatrix[,3],
  cat.mu = covMatrix[,4],
  cat.sd = covMatrix[,5],
  fox.mu = covMatrix[,8],
  fox.sd = covMatrix[,9],
  nNight = max(n_long$night),
  moon = moonArray,
  jDate = jDateArray,
  effort = effortArray,
  y = my_array, 
  nSP = max(n_long$Season),
  season = seasonData,
  nseason = 3
)
save(data_list, file = "./results/data_list.Rdata") #not really a 'result' but saving there for later use

# initial values for N should be maximum possible counts
N_init = array(10, dim = c(data_list$nsite, data_list$nseason))
# set to NA for all seasons t>1 because we cannot know this (needs to be specified by R + S)
N_init[,2:6] <- NA
# initial values for recruits can also be set to maximum possible counts
R_init = array(10, dim = c(data_list$nsite, data_list$nseason))
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
system("say Calculations Complete.")
saveRDS(my_mod, "./results/rod_model.RDS")
summary(my_mod, vars = c("beta0","beta1",
                         "phi0", "phi1", "phi2", "phi3", "phi4", "phi5", "phi6",
                         "gamma0", "gamma1", "gamma2",
                         "alpha0", "alpha1", "alpha2", "alpha3"
))
plot(my_mod)


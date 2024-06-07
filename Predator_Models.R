################################################################################
#                                                                              #
#    Multi-Species Autologistic Occupancy Model for Determining Average        #
#    Site Occupancy of Carnivores in Iowa City, IA                             #
#                                                                              #
#    Last Updated: 07 Jun 2024                                                 #
#                                                                              #
################################################################################
# This script prepares the data and runs an autologistic multi-species occupancy
# model for generating estimates of average site occupancy probability for
# domestic cats, coyotes, red fox, and American mink to use in the model of
# deer mouse abundance

#### Load & Clean Data for Model ####
# load functions used to clean data
functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}
# read in the full predator dataset with all species
dat <- read.csv("./data/predatorOccu.csv")
# converting the site names & species into unique numbers
dat$Site <- as.numeric(as.factor(dat$Site))
dat$Species <- as.numeric(as.factor(dat$Species))
str(dat)
# converting Y to NA if J (i.e., trap-days) is 0
dat$Y[dat$J==0] <- NA
# re-order the data to be in order by species, season, then site
dat <- dat[order(
  dat[,"Species"], 
  dat[,"Season"],
  dat[,"Site"]), ]
# create blank array w/ correct dimensions then fill it in w/ detections
y <- array(data = NA, dim = c(max(dat$Species), max(dat$Site), max(dat$Season)))
for(i in 1:nrow(dat)){
  y[
    dat$Species[i],
    dat$Site[i],
    dat$Season[i]
  ] <- dat$Y[i]
}
# create blank matrix w/ correct dimensions then fill it in w/ camera-Days
J <- matrix(data = NA, nrow = max(dat$Site), ncol = max(dat$Season))
for(i in 1:nrow(dat)){
  J[
    dat$Site[i],
    dat$Season[i]
  ] <- dat$J[i]
}

# creating a dummy 'season' variable for the detection submodel
seasonData <- as.factor(c("fall","winter","spring","summer",
                          "fall","winter","spring","summer",
                          "fall","winter","spring","summer",
                          "fall"))

# reading in site covariates & converting sites to numeric factors
siteCovs <- read.csv("./data/siteCovs_Pred.csv")
siteCovs$Site <- as.numeric(as.factor(siteCovs$Site))
cor(scale(siteCovs))
# re-order the data to be in order by site
siteCovs <- siteCovs[order(siteCovs[,"Site"]), ]
siteCovs$Imperv <- as.numeric(scale(siteCovs$Imperv))
siteCovs$Forest <- as.numeric(scale(siteCovs$Forest))
siteCovs$Prairie <- as.numeric(scale(siteCovs$Prairie))
siteCovs$ResUnits <- as.numeric(scale(siteCovs$ResUnits))
siteCovs$Crop <- as.numeric(scale(siteCovs$Crop))
siteCovs$Dist_to_Wat <- as.numeric(scale(siteCovs$Dist_to_Wat))

# last but not least, making the array of possible occupancy states
# (e.g., 1/1/1/1 = all species present; 0/1/0/0 = only coyote present, etc.)
Xcat <- matrix(c(1, 1, 1, 1,
                 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1,
                 1, 1, 0, 0,
                 1, 0, 1, 0,
                 1, 0, 0, 1,
                 0, 1, 1, 0,
                 0, 1, 0, 1,
                 0, 0, 1, 1,
                 0, 0, 0, 0), ncol = 4, byrow = TRUE)

#### RUN MODEL ####
# Data list for model
data_list <- list(
  nspec = max(dat$Species),
  nsite = max(dat$Site),
  house = siteCovs$ResUnits,
  crop = siteCovs$Crop,
  forest = siteCovs$Forest,
  prairie = siteCovs$Prairie,
  water = siteCovs$Dist_to_Wat,
  imperv = siteCovs$Imperv,
  season = seasonData,
  nseason = 4,
  y = y, 
  J = J,
  Xcat = Xcat,
  nSP = max(dat$Season)
)
# Fit the model
library(runjags)
pred_mod <- runjags::run.jags(
  model = "./jags/predatorModel.R",
  monitor = c("a0", "a1", "a2", "b0", "b1", "b2", "c0", "c1", "d0", "d1", "d2",
              "e0", "e1", "g0", "g1", "h0", "h1",
              "l0", "l1", "m0", "m1",
              "n0", "n1",
              "p0", "p1",
              "phi",
              "x"),
  data = data_list,
  n.chains = 3,
  inits = pred_inits,
  burnin = 1000,
  sample = 5000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
# Save output for later use & print model summary to check convergence, etc.
saveRDS(pred_mod, "./results/pred_model.RDS")
summary(pred_mod,
        vars = c("a0", "a1", "a2", "b0", "b1", "b2", "c0", "c1", "d0", "d1", "d2",
                "e0", "e1", "g0", "g1", "h0", "h1", "l0", "l1", "m0", "m1", "n0", "n1",
                "p0", "p1",
                "phi"
        ))
plot(pred_mod)

#### Post-Model Calculations ####
# convert the model file to an mcmc object; use MCMCvis to get summary statistics
library(coda)
mc <- as.mcmc(pred_mod)
library(MCMCvis)
MCMCvis::MCMCsummary(mc,
                     params = c("a0", "a1", "a2", "b0", "b1", "b2", "c0", "c1", "d0", "d1", "d2",
                                "e0", "e1", "g0", "g1", "h0", "h1", "l0", "l1", "m0", "m1", "n0", "n1",
                                "p0", "p1",
                                "phi"
                                ),
                     probs = c(0.025, 0.5, 0.975), 
                     round = 2)

# convert the mcmc object to a matrix
mc <- as.matrix(mc)
# sub-sample the mcmc matrix a bit as we don't really need to make predictions with 
# all 60K samples
set.seed(554)
mc_sub <- mc[sample(1:nrow(mc), 10000), ]
# and use split_mcmc
mc <- split_mcmc(mc_sub)
rm(mc_sub)

# the script below uses the model output to predict which occupancy state (e.g., 
# cats and coyotes present, only fox present, etc) of each site at the 1st timestep.
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(1),
    12
  )
)
# no species
pred_psi[,1] <- 1
# Species 1 (cat)
pred_psi[,2] <- exp(mc$a0[,1])
# Species 2 (coyote)
pred_psi[,3] <- exp(mc$b0[,1])
# Species 3 (fox)
pred_psi[,4] <- exp(mc$c0[,1])
# Species 4 (mink)
pred_psi[,5] <- exp(mc$d0[,1])
# species 1 & 2
pred_psi[,6] <- exp(mc$a0[,1] + mc$b0[,1] + mc$e0[,1])
# species 1 & 3
pred_psi[,7] <- exp(mc$a0[,1] + mc$c0[,1] + mc$g0[,1])
# species 1 & 4
pred_psi[,8] <- exp(mc$a0[,1] + mc$d0[,1] + mc$h0[,1])
# species 2 & 3
pred_psi[,9] <- exp(mc$b0[,1] + mc$c0[,1] + mc$l0[,1])
# species 2 & 4
pred_psi[,10] <- exp(mc$b0[,1] + mc$d0[,1] + mc$m0[,1])
# species 3 & 4
pred_psi[,11] <- exp(mc$c0[,1] + mc$d0[,1] + mc$n0[,1])
# all species together
pred_psi[,12] <- exp(mc$a0[,1] + mc$b0[,1] + mc$c0[,1] + mc$d0[,1] + 
                       mc$e0[,1] + mc$g0[,1] + mc$h0[,1] + mc$l0[,1] + mc$m0[,1] + mc$n0[,1])
# convert to probability
prob_psi <- t(apply(pred_psi, 1, function(x) x/sum(x)))
# the following code does the same thing as the above, but with the addition of the autologistic
# term to account for temporal correlation (thus calculating occupancy state in the subsequent
# timesteps
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(1),
    12
  )
)
# no species
pred_psi_auto[,1] <- 1
pred_psi_auto[,2] <- exp(mc$a0[,1] + mc$phi[,1])
pred_psi_auto[,3] <- exp(mc$b0[,1] + mc$phi[,2])
pred_psi_auto[,4] <- exp(mc$c0[,1] + mc$phi[,3])
pred_psi_auto[,5] <- exp(mc$d0[,1] + mc$phi[,4])
pred_psi_auto[,6] <- exp(mc$a0[,1] + mc$phi[,1] + mc$b0[,1] + mc$phi[,2] + mc$e0[,1])
pred_psi_auto[,7] <- exp(mc$a0[,1] + mc$phi[,1] + mc$c0[,1] + mc$phi[,3] + mc$g0[,1])
pred_psi_auto[,8] <- exp(mc$a0[,1] + mc$phi[,1] + mc$d0[,1] + mc$phi[,4] + mc$h0[,1])
pred_psi_auto[,9] <- exp(mc$b0[,1] + mc$phi[,2] + mc$c0[,1] + mc$phi[,3] + mc$l0[,1])
pred_psi_auto[,10] <- exp(mc$b0[,1] + mc$phi[,2] + mc$d0[,1] + mc$phi[,4] + mc$m0[,1])
pred_psi_auto[,11] <- exp(mc$c0[,1] + mc$phi[,3] + mc$d0[,1] + mc$phi[,4] + mc$n0[,1])
pred_psi_auto[,12] <- exp(mc$a0[,1] + mc$phi[,1] + mc$b0[,1] + mc$phi[,2] + mc$c0[,1] + 
                            mc$phi[,3] + mc$d0[,1] + mc$phi[,4] +
                       mc$e0[,1] + mc$g0[,1] + mc$h0[,1] + mc$l0[,1] + mc$m0[,1] + mc$n0[,1])
# convert to probability
prob_psi_auto <- t(apply(pred_psi, 1, function(x) x/sum(x)))

# calculating average psi (i.e., average occupancy probabilty across time)
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)
trueProb <- apply(
  trueProb,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Marginal Occupancy ####
# calculates the marginal occupancy probability (i.e., irrespective of the occupancy of any other species) for each
# species in the model
Cat <- trueProb[,2] + trueProb[,6] + trueProb[,7] + trueProb[,8] + trueProb[,12]
Coyote <- trueProb[,3] + trueProb[,6] + trueProb[,9] + trueProb[,10] + trueProb[,12]
Fox <- trueProb[,4] + trueProb[,7] + trueProb[,9] + trueProb[,11] + trueProb[,12]
Mink <- trueProb[,5] + trueProb[,8] + trueProb[,10] + trueProb[,11] + trueProb[,12]
#### Conditional Occupancies ####
# calculate the probability of the small-case species given the presence/absence of the Big Case species
catCoyote <- (trueProb[,6] + trueProb[,12]) / (trueProb[,6] + trueProb[,12] + trueProb[,3] + trueProb[,9] + trueProb[,10])
# then we compare modeled co-occurrence (smallBig) with co-occurrence if the species occurred together randomly (small*Big)
catCoyote
Cat*Coyote # cats & coyotes occur together MORE than expected
# more likely to co-occur as imperv increases
# likely b/c of increased probability of occurrence of cats?
foxCoyote <- (trueProb[,9] + trueProb[,12]) / (trueProb[,9] + trueProb[,12] + trueProb[,3] + trueProb[,6] + trueProb[,10])
foxCoyote
Fox*Coyote # fox & coyotes occur together MORE than expected
# less likely to co-occur as imperv increases
# likely b/c foxes are more likely to occur as imperv increases, while coyotes avoid imperv
minkCoyote <- (trueProb[,10] + trueProb[,12]) / (trueProb[,10] + trueProb[,12] + trueProb[,3] + trueProb[,6] + trueProb[,9])
minkCoyote
Mink*Coyote # mink & coyotes occur together AS EXPECTED

#### Average Occupancy by Site ####
# we're going to calculate the average occupancy probability through time at each site for each species
# first we generate a blank array to later hold the predictions
x_array <- as.array(mc$x)
pred_Site <- array(
  NA,
  dim = c(
    2,
    max(data_list$nsite),
    max(data_list$nspec)
  )
)
# then we take the predicted occupancy states (1,2,3,4,etc.) and convert them to a '0' if the state
# does not include our species of interest (e.g., for cats, occupancy states 3 [fox only] and
# 4 [mink only] do not include cats) and a '1' if it does
# cat marginal occupancy
cat <- replace(x_array, x_array == 12, 0)
cat <- replace(cat, cat == 11, 0)
cat <- replace(cat, cat == 10, 0)
cat <- replace(cat, cat == 9, 0)
cat <- replace(cat, cat == 5, 0)
cat <- replace(cat, cat == 4, 0)
cat <- replace(cat, cat == 3, 0)
cat <- replace(cat, cat == 2, 1)
cat <- replace(cat, cat == 6, 1)
cat <- replace(cat, cat == 7, 1)
cat <- replace(cat, cat == 8, 1)
pred_Site[1,,1] <- apply(cat, 2, mean)
pred_Site[2,,1] <- apply(cat, 2, sd)
# coyote marginal occupancy
coyote <- replace(x_array, x_array == 12, 0)
coyote <- replace(coyote, coyote == 11, 0)
coyote <- replace(coyote, coyote == 8, 0)
coyote <- replace(coyote, coyote == 7, 0)
coyote <- replace(coyote, coyote == 5, 0)
coyote <- replace(coyote, coyote == 4, 0)
coyote <- replace(coyote, coyote == 2, 0)
coyote <- replace(coyote, coyote == 3, 1)
coyote <- replace(coyote, coyote == 6, 1)
coyote <- replace(coyote, coyote == 9, 1)
coyote <- replace(coyote, coyote == 10, 1)
pred_Site[1,,2] <- apply(coyote, 2, mean)
pred_Site[2,,2] <- apply(coyote, 2, sd)
# fox marginal occupancy
fox <- replace(x_array, x_array == 12, 0)
fox <- replace(fox, fox == 10, 0)
fox <- replace(fox, fox == 8, 0)
fox <- replace(fox, fox == 6, 0)
fox <- replace(fox, fox == 5, 0)
fox <- replace(fox, fox == 2, 0)
fox <- replace(fox, fox == 3, 0)
fox <- replace(fox, fox == 4, 1)
fox <- replace(fox, fox == 7, 1)
fox <- replace(fox, fox == 9, 1)
fox <- replace(fox, fox == 11, 1)
pred_Site[1,,3] <- apply(fox, 2, mean)
pred_Site[2,,3] <- apply(fox, 2, sd)
# mink marginal occupancy
mink <- replace(x_array, x_array == 12, 0)
mink <- replace(mink, mink == 9, 0)
mink <- replace(mink, mink == 7, 0)
mink <- replace(mink, mink == 6, 0)
mink <- replace(mink, mink == 2, 0)
mink <- replace(mink, mink == 3, 0)
mink <- replace(mink, mink == 4, 0)
mink <- replace(mink, mink == 5, 1)
mink <- replace(mink, mink == 8, 1)
mink <- replace(mink, mink == 10, 1)
mink <- replace(mink, mink == 11, 1)
pred_Site[1,,4] <- apply(mink, 2, mean)
pred_Site[2,,4] <- apply(mink, 2, sd)
# transpose dataframe
c <- pred_Site[,,1]
c <- t(c)
y <- pred_Site[,,2]
y <- t(y)
x <- pred_Site[,,3]
x <- t(x)
m <- pred_Site[,,4]
m <- t(m)
# merge all of them together and save as a dataframe for later use
occuSitePreds <- cbind(siteCovs$Site,c,y,x,m)
colnames(occuSitePreds) <- c("camera","cat.mu","cat.sd","coyote.mu","coyote.sd",
                             "fox.mu","fox.sd","mink.mu","mink.sd")
write.csv(occuSitePreds, "./data/modeledOccu_Predators.csv")

#### Additional Calculations (Cat Detection by Season, Coyote/Fox Co-occurrence & Impervious Relationship ####
# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
predV <- c(1,2,3,4)
pred_mat <- cbind(
  1,
  predV
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# create the matrix of the coefficients from the model
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$p0),
    ncol(pred_mat)
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_p <- cbind(mc$p0[,2],mc$p1[,2]) %*% pred_mat

# convert to probability
prob_p <- apply(pred_p,
                c(2),
                logit2prob)
prob_p <- apply(
  prob_p,
  c(2),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

#### coyote + fox imperv relationship ####
predV <- seq(-1.15, 2.47,length.out=200)
pred_mat <- cbind(
  1,
  predV
)
pred_mat <- t(pred_mat)

pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    12
  )
)

# no species
pred_psi[,,1] <- 1
# Species 1 (cat)
pred_psi[,,2] <- exp(mc$a0[,1])
# Species 2 (coyote)
pred_psi[,,3] <- exp(mc$b0[,1])
# Species 3 (fox)
pred_psi[,,4] <- exp(mc$c0[,1])
# Species 4 (mink)
pred_psi[,,5] <- exp(mc$d0[,1])
# species 1 & 2
pred_psi[,,6] <- exp(mc$a0[,1] + mc$b0[,1] + 
                       cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat)
# species 1 & 3
pred_psi[,,7] <- exp(mc$a0[,1] + mc$c0[,1] + 
                       cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat)
# species 1 & 4
pred_psi[,,8] <- exp(mc$a0[,1] + mc$d0[,1] + 
                       cbind(mc$h0[,1],mc$h1[,1]) %*% pred_mat)
# species 2 & 3
pred_psi[,,9] <- exp(mc$b0[,1] + mc$c0[,1] + 
                      cbind(mc$l0[,1],mc$l0[,1]) %*% pred_mat)
# species 2 & 4
pred_psi[,,10] <- exp(mc$b0[,1] + mc$d0[,1] + 
                       cbind(mc$m0[,1],mc$m1[,1]) %*% pred_mat)
# species 3 & 4
pred_psi[,,11] <- exp(mc$c0[,1] + mc$d0[,1] + 
                       cbind(mc$n0[,1],mc$n1[,1]) %*% pred_mat)
# all species together
pred_psi[,,12] <- exp(mc$a0[,1] + mc$b0[,1] + mc$c0[,1] + mc$d0[,1] + 
                        cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat + 
                        cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat + 
                        cbind(mc$h0[,1],mc$h1[,1]) %*% pred_mat + 
                        cbind(mc$l0[,1],mc$l0[,1]) %*% pred_mat + 
                        cbind(mc$m0[,1],mc$m1[,1]) %*% pred_mat + 
                        cbind(mc$n0[,1],mc$n1[,1]) %*% pred_mat)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# no species
pred_psi_auto[,,1] <- 1
pred_psi_auto[,,2] <- exp(mc$a0[,1] + mc$phi[,1])
pred_psi_auto[,,3] <- exp(mc$b0[,1] + mc$phi[,2])
pred_psi_auto[,,4] <- exp(mc$c0[,1] + mc$phi[,3])
pred_psi_auto[,,5] <- exp(mc$d0[,1] + mc$phi[,4])
pred_psi_auto[,,6] <- exp(mc$a0[,1] + mc$phi[,1] + mc$b0[,1] + mc$phi[,2] + 
                            cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat)
pred_psi_auto[,,7] <- exp(mc$a0[,1] + mc$phi[,1] + mc$c0[,1] + mc$phi[,3] + 
                            cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat)
pred_psi_auto[,,8] <- exp(mc$a0[,1] + mc$phi[,1] + mc$d0[,1] + mc$phi[,4] + 
                            cbind(mc$h0[,1],mc$h1[,1]) %*% pred_mat)
pred_psi_auto[,,9] <- exp(mc$b0[,1] + mc$phi[,2] + mc$c0[,1] + mc$phi[,3] + 
                            cbind(mc$l0[,1],mc$l1[,1]) %*% pred_mat)
pred_psi_auto[,,10] <- exp(mc$b0[,1] + mc$phi[,2] + mc$d0[,1] + mc$phi[,4] + 
                             cbind(mc$m0[,1],mc$m1[,1]) %*% pred_mat)
pred_psi_auto[,,11] <- exp(mc$c0[,1] + mc$phi[,3] + mc$d0[,1] + mc$phi[,4] + 
                             cbind(mc$n0[,1],mc$n1[,1]) %*% pred_mat)
pred_psi_auto[,,12] <- exp(mc$a0[,1] + mc$phi[,1] + mc$b0[,1] + mc$phi[,2] + mc$c0[,1] + 
                            mc$phi[,3] + mc$d0[,1] + mc$phi[,4] +
                             cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat + 
                             cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat + 
                             cbind(mc$h0[,1],mc$h1[,1]) %*% pred_mat + 
                             cbind(mc$l0[,1],mc$l1[,1]) %*% pred_mat + 
                             cbind(mc$m0[,1],mc$m1[,1]) %*% pred_mat + 
                             cbind(mc$n0[,1],mc$n1[,1]) %*% pred_mat)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

trueProb <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

coyFox <- rbind(data.frame(imperv = seq(-1.15, 2.47,length.out=200),
                            psi = ((trueProb[3,,3] + trueProb[3,,6] + trueProb[3,,9] + trueProb[3,,10] + trueProb[3,,12])*(trueProb[3,,4] + trueProb[3,,7] + trueProb[3,,9] + trueProb[3,,11] + trueProb[3,,12])),
                            upper = ((trueProb[5,,3] + trueProb[5,,6] + trueProb[5,,9] + trueProb[5,,10] + trueProb[5,,12])*(trueProb[5,,4] + trueProb[5,,7] + trueProb[5,,9] + trueProb[5,,11] + trueProb[5,,12])),
                            lower = ((trueProb[1,,3] + trueProb[1,,6] + trueProb[1,,9] + trueProb[1,,10] + trueProb[1,,12])*(trueProb[1,,4] + trueProb[1,,7] + trueProb[1,,9] + trueProb[1,,11] + trueProb[1,,12])),
                            state = "Expected Co-occurrence"),
                 data.frame(imperv = seq(-1.15, 2.47,length.out=200),
                            psi = (trueProb[3,,9] + trueProb[3,,12]) / (trueProb[3,,3] + trueProb[3,,6] + trueProb[3,,9] + trueProb[3,,10] + trueProb[3,,12]),
                            lower = (trueProb[1,,9] + trueProb[1,,12]) / (trueProb[1,,3] + trueProb[1,,6] + trueProb[1,,9] + trueProb[1,,10] + trueProb[1,,12]),
                            upper = (trueProb[5,,9] + trueProb[5,,12]) / (trueProb[5,,3] + trueProb[5,,6] + trueProb[5,,9] + trueProb[5,,10] + trueProb[5,,12]),
                            state = "Actual Co-occurrence")
)

ggplot(coyFox, aes(x=imperv, y=psi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=state), alpha=.3, show.legend=FALSE) +
  scale_fill_manual(values=c("gray30","gray60"), guide = "none") +
  geom_line(aes(color=state, linetype=state), linewidth=1.05, show.legend = FALSE) +
  scale_color_manual(values=c("gray30","gray60"), guide = "none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(title="", x="", y="", linetype="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=10))

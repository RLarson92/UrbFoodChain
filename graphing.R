library(coda)
library(runjags)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(grafify)

# load functions used to clean data
functions_to_load <- list.files(
  "./functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}

my_mod <- readRDS("./results/rod_model.RDS")
load("./results/data_list.Rdata")

mc <- as.mcmc(my_mod, chains = FALSE) # need both 'runjags' and 'coda' loaded for this to work
# convert the mcmc object to a matrix
mc <- as.matrix(mc)
# sub-sample the mcmc matrix a bit as we don't really need to make predictions with 
# all 60K samples
set.seed(554)
mc_sub <- mc[sample(1:nrow(mc), 10000), ]
# and use split_mcmc
mc <- split_mcmc(mc_sub)
rm(mc_sub)

range(data_list$PC1)
mean(data_list$contag)
##### Persistence (phi[j,t]) #####
##### Cats #####
# grab quantiles of cat occupancy across all sites
quantile(mc$cat)
# 0%          25%          50%          75%         100% 
# 0.00       0.08         0.21         0.35          0.99
predV_cat3 <- cbind(1,                                              # intercept                     
                      seq(-3.77, 1.92,length.out=200),              # PC1
                      1,                                            # cats present
                      mean(mc$fox),                                 # fox; set to mean occupancy across all sites
                      seq(-3.77, 1.92,length.out=200)*0.75,         # PC1*cat interaction
                      seq(-3.77, 1.92,length.out=200)*mean(mc$fox), # PC1*fox interaction
                    mean(data_list$contag)                          # contagion; set to mean
                      )
predV_cat3 <- t(predV_cat3)
# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$phi0),
    ncol(predV_cat3)
  )
)
pred_phi <- cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1],mc$phi3[,1],mc$phi4[,1],mc$phi5[,1],mc$phi6[,1]) %*% predV_cat3
prob_phi <- apply(pred_phi,
                  c(2),
                  logit2prob)
truePhi_cat3 <- apply(
  prob_phi,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

predV_cat2 <- cbind(1,                                              # intercept                     
                    seq(-3.77, 1.92,length.out=200),              # PC1
                    1,                                            # cats present
                    mean(mc$fox),                                 # fox; set to mean occupancy across all sites
                    seq(-3.77, 1.92,length.out=200)*0.5,         # PC1*cat interaction
                    seq(-3.77, 1.92,length.out=200)*mean(mc$fox), # PC1*fox interaction
                    mean(data_list$contag)                          # contagion; set to mean
)
predV_cat2 <- t(predV_cat2)
# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$phi0),
    ncol(predV_cat2)
  )
)
pred_phi <- cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1],mc$phi3[,1],mc$phi4[,1],mc$phi5[,1],mc$phi6[,1]) %*% predV_cat2
prob_phi <- apply(pred_phi,
                  c(2),
                  logit2prob)
truePhi_cat2 <- apply(
  prob_phi,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

predV_cat1 <- cbind(1,                                            # intercept                     
                    seq(-3.77, 1.92,length.out=200),              # PC1
                    0,                                            # cats absent
                    mean(mc$fox),                                 # fox; set to mean occupancy across all sites
                    seq(-3.77, 1.92,length.out=200)*0.25,         # PC1*cat interaction
                    seq(-3.77, 1.92,length.out=200)*mean(mc$fox), # PC1*fox interaction
                    mean(data_list$contag)                        # contagion; set to mean
)
predV_cat1 <- t(predV_cat1)
# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$phi0),
    ncol(predV_cat1)
  )
)
pred_phi <- cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1],mc$phi3[,1],mc$phi4[,1],mc$phi5[,1],mc$phi6[,1]) %*% predV_cat1
prob_phi <- apply(pred_phi,
                  c(2),
                  logit2prob)
truePhi_cat1 <- apply(
  prob_phi,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

datC <- rbind(data.frame(group = "Pr = 0.75",
                         PC1 = seq(-3.77, 1.92,length.out=200),
                         phi = truePhi_cat3[3,],
                         upper = truePhi_cat3[5,],
                         lower = truePhi_cat3[1,]),
              data.frame(group = "Pr = 0.50",
                         PC1 = seq(-3.77, 1.92,length.out=200),
                         phi = truePhi_cat2[3,],
                         upper = truePhi_cat2[5,],
                         lower = truePhi_cat2[1,]),
              data.frame(group = "Pr = 0.25",
                         PC1 = seq(-3.77, 1.92,length.out=200),
                         phi = truePhi_cat0[3,],
                         upper = truePhi_cat0[5,],
                         lower = truePhi_cat0[1,])
              )
p1 <- ggplot(datC, aes(x=PC1, y=phi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.3, show.legend=FALSE) +
  scale_fill_manual(values = c("#5cae63","#999999","#986eac")) +
  geom_line(aes(color=group), show.legend=FALSE) +
  scale_color_manual(values = c("#5cae63","#999999","#986eac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="", y="Persistence Probability", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title.y=element_text(size=9), axis.title.x = element_blank())

##### Foxes #####
quantile(mc$fox)
# 0%         25%         50%         75%        100% 
# 0.00      0.15        0.27        0.46        0.99 
predV_fox3 <- cbind(1,                                            # intercept                     
                    seq(-3.77, 1.92,length.out=200),              # PC1
                    mean(mc$cat),                                 # cat; set to mean occupancy across all sites
                    1,                                            # fox; present
                    seq(-3.77, 1.92,length.out=200)*mean(mc$cat), # PC1*cat interaction
                    seq(-3.77, 1.92,length.out=200)*0.75,            # PC1*fox interaction
                    mean(data_list$contag)                        # contagion; set to mean
)
predV_fox3 <- t(predV_fox3)
# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$phi0),
    ncol(predV_fox3)
  )
)
pred_phi <- cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1],mc$phi3[,1],mc$phi4[,1],mc$phi5[,1],mc$phi6[1]) %*% predV_fox3
prob_phi <- apply(pred_phi,
                  c(2),
                  logit2prob)
truePhi_fox3 <- apply(
  prob_phi,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

predV_fox2 <- cbind(1,                                            # intercept                     
                    seq(-3.77, 1.92,length.out=200),              # PC1
                    mean(mc$cat),                                 # cat; set to mean occupancy across all sites
                    0,                                            # fox; absent
                    seq(-3.77, 1.92,length.out=200)*mean(mc$cat), # PC1*cat interaction
                    seq(-3.77, 1.92,length.out=200)*0.5,            # PC1*fox interaction
                    mean(data_list$contag)                        # contagion; set to mean
)
predV_fox2 <- t(predV_fox2)
# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$phi0),
    ncol(predV_fox2)
  )
)
pred_phi <- cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1],mc$phi3[,1],mc$phi4[,1],mc$phi5[,1],mc$phi6[,1]) %*% predV_fox2
prob_phi <- apply(pred_phi,
                  c(2),
                  logit2prob)
truePhi_fox2 <- apply(
  prob_phi,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

predV_fox1 <- cbind(1,                                            # intercept                     
                    seq(-3.77, 1.92,length.out=200),              # PC1
                    mean(mc$cat),                                 # cat; set to mean occupancy across all sites
                    0,                                            # fox; absent
                    seq(-3.77, 1.92,length.out=200)*mean(mc$cat), # PC1*cat interaction
                    seq(-3.77, 1.92,length.out=200)*0.25,            # PC1*fox interaction
                    mean(data_list$contag)                        # contagion; set to mean
)
predV_fox1 <- t(predV_fox1)
# create the matrix of the coefficients from the model
pred_phi <- array(
  NA,
  dim = c(
    nrow(mc$phi0),
    ncol(predV_fox1)
  )
)
pred_phi <- cbind(mc$phi0[,1],mc$phi1[,1],mc$phi2[,1],mc$phi3[,1],mc$phi4[,1],mc$phi5[,1],mc$phi6[,1]) %*% predV_fox1
prob_phi <- apply(pred_phi,
                  c(2),
                  logit2prob)
truePhi_fox1 <- apply(
  prob_phi,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

datF <- rbind(data.frame(group = "Pr = 0.75",
                        PC1 = seq(-3.77, 1.92,length.out=200),
                        phi = truePhi_fox3[3,],
                        upper = truePhi_fox3[5,],
                        lower = truePhi_fox3[1,]),
              data.frame(group = "Pr = 0.50",
                         PC1 = seq(-3.77, 1.92,length.out=200),
                         phi = truePhi_fox2[3,],
                         upper = truePhi_fox2[5,],
                         lower = truePhi_fox2[1,]),
             data.frame(group = "Pr = 0.25",
                        PC1 = seq(-3.77, 1.92,length.out=200),
                        phi = truePhi_fox1[3,],
                        upper = truePhi_fox1[5,],
                        lower = truePhi_fox1[1,])
)
p2<-ggplot(datF, aes(x=PC1, y=phi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.3, show.legend=FALSE) +
  scale_fill_manual(values = c("#5cae63","#999999","#986eac")) +
  geom_line(aes(color=group), show.legend=FALSE) +
  scale_color_manual(values = c("#5cae63","#999999","#986eac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="", y="", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title=element_blank())

# empty plot for legend
p4 <- ggplot(datF, aes(x=PC1, y=phi))+
  geom_line(aes(color=group)) +
  scale_color_manual(values = c("#5cae63","#999999","#986eac"), name = "Predator\nOccupancy\nProbability") +
  theme_bw() +
  theme(legend.position = c(0.5,0.5), legend.title = element_text(size=9), legend.background=element_blank(),
        legend.text = element_text(size=6))
leg <- ggpubr::get_legend(p4)
leg <- ggpubr::as_ggplot(leg)

jpeg("./results/phi.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2,leg, ncol=3, widths=c(2,1.89,0.75),
                        bottom=textGrob("PC1", gp=gpar(fontsize=9)))
dev.off()

##### Recruitment (gamma[j,t]) #####
range(data_list$contag)
# spring
predV_SP <- cbind(1,                                # intercept                     
               seq(-1.3, 2.164393,length.out=200),  # contagion
               1                                    # season; 2 = spring
)
predV_SP <- t(predV_SP)
# create the matrix of the coefficients from the model
pred_gamma <- array(
  NA,
  dim = c(
    nrow(mc$gamma0),
    ncol(predV_SP)
  )
)
pred_gamma <- cbind(mc$gamma0[,1],mc$gamma1[,1],mc$gamma2[,2]) %*% predV_SP
prob_gamma <- apply(pred_gamma,
                  c(2),
                  exp)
trueGammaSP <- apply(
  prob_gamma,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

# summer
predV_SU <- cbind(1,                                   # intercept                     
                  seq(-1.3, 2.164393,length.out=200),  # contagion
                  1                                    # season; 3 = summer
)
predV_SU <- t(predV_SU)
# create the matrix of the coefficients from the model
pred_gamma <- array(
  NA,
  dim = c(
    nrow(mc$gamma0),
    ncol(predV_SU)
  )
)
pred_gamma <- cbind(mc$gamma0[,1],mc$gamma1[,1],mc$gamma2[,3]) %*% predV_SU
prob_gamma <- apply(pred_gamma,
                    c(2),
                    exp)
trueGammaSU <- apply(
  prob_gamma,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

# fall
predV_FA <- cbind(1,                                   # intercept                     
                  seq(-1.3, 2.164393,length.out=200),  # contagion
                  1                                    # season; 1 = fall
)
predV_FA <- t(predV_FA)
# create the matrix of the coefficients from the model
pred_gamma <- array(
  NA,
  dim = c(
    nrow(mc$gamma0),
    ncol(predV_FA)
  )
)
pred_gamma <- cbind(mc$gamma0[,1],mc$gamma1[,1],mc$gamma2[,1]) %*% predV_FA
prob_gamma <- apply(pred_gamma,
                    c(2),
                    exp)
trueGammaFA <- apply(
  prob_gamma,
  2,
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975),
  (na.rm = TRUE)
)

dat <- rbind(data.frame(group = "Spring",
                         contag = seq(13, 100,length.out=200),
                         gamma = trueGammaSP[3,],
                         upper = trueGammaSP[5,],
                         lower = trueGammaSP[1,]),
              data.frame(group = "Summer",
                         contag = seq(13, 100,length.out=200),
                         gamma = trueGammaSU[3,],
                         upper = trueGammaSU[5,],
                         lower = trueGammaSU[1,]),
             data.frame(group = "Fall",
                        contag = seq(13, 100,length.out=200),
                        gamma = trueGammaFA[3,],
                        upper = trueGammaFA[5,],
                        lower = trueGammaFA[1,])
)
jpeg("./results/gamma.jpeg", width = 6, height = 4, units = 'in', res = 300)
ggplot(dat, aes(x=contag, y=gamma)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.3, show.legend=FALSE) +
  scale_fill_manual(values = c("#ccbb44","#ee6677","#228833")) +
  geom_line(aes(color=group)) +
  scale_color_manual(values = c("#ccbb44","#ee6677","#228833")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  labs(x="Contagion", y="Recruitment", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.position=c(0.08,0.95), legend.background=element_blank()) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
dev.off()

##### Capture Rate (Supplemental Figure) #####
##### Moon Illumination #####
# generate a sequence of moon illumination values
predV <- seq(-1.343997, 1.474745, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  predV,
  0,  # mean of scale(jDate)
  0   # mean of scale(effort)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# create the matrix of the coefficients from the model
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat)
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_p <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

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

## Graphing
moon <- rbind(data.frame(moon = seq(0, 1, length.out = 200),
                         p = prob_p[3,],
                         upper = prob_p[5,],
                         lower = prob_p[1,])
)

library(scales)
p1 <- ggplot(moon, aes(x=moon, y=p)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.5, show.legend=FALSE) +
  geom_line(aes(), show.legend=FALSE) +
  scale_x_continuous(labels = label_number(accuracy = 0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="Moon Illumination\n(proportion full)", y="Capture Probability", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9))
##### Julian Date #####
# generate a sequence of Julian date values
# [1] -1.551832  1.596653
predV <- seq(-1.6, 1.6, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  0,  # mean of scale(moon)
  predV,
  0   # mean of scale(effort)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# create the matrix of the coefficients from the model
pred_p <- array(
  NA,
  dim = c(
    nrow(mc$alpha0),
    ncol(pred_mat)
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_p <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

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

## Graphing
range(jDate_long$jDate)
# [1] 101 310

jDate <- rbind(data.frame(jDate = seq(101, 310, length.out = 200),
                          p = prob_p[3,],
                          upper = prob_p[5,],
                          lower = prob_p[1,])
)

p2 <- ggplot(jDate, aes(x=jDate, y=p)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.5, show.legend=FALSE) +
  geom_line(aes(), show.legend=FALSE) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="Ordinal Date\n", y="", color="", linewidth="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9))
##### Effort #####
# generate a sequence of moon illumination values
# [1] -2.544777  1.402502
predV <- c(-2.5447767, -2.1500488, -1.7553209, -1.3605930, -0.9658651, 
           -0.5711372, -0.1764093, 0.2183186, 0.6130466, 1.0077745, 1.4025024)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# contag values, and adding a row of 1's for the intercept & other predictor values
pred_mat <- cbind(
  1,
  0,     # mean of scale(moon)
  0,     # mean of scale(temp)
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
    nrow(mc$alpha0),
    ncol(pred_mat)
  )
)

# dimensions should be: number of MCMC interations by length of sequence of predictor variable
# by number of species
# fill it in
pred_p <- cbind(mc$alpha0[,1],mc$alpha1[,1],mc$alpha2[,1],mc$alpha3[,1]) %*% pred_mat

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

## Graphing
effort <- rbind(data.frame(effort = c(0,1,2,3,4,5,6,7,8,9,10),
                           p = prob_p[3,],
                           upper = prob_p[5,],
                           lower = prob_p[1,])
                )

p3<-ggplot(effort, aes(x=effort, y=p)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), show.legend=FALSE) +
  geom_point() +
  scale_color_manual(values = c("#000000")) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  labs(x="Trap Effort\n", y="", size="", color="Group") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.title=element_blank(), legend.position=c(0.6,0.875), legend.background=element_blank(),
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))

jpeg("./results/p.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2,p3, ncol=3)
dev.off()

##### Supplemental Figure #####
cat <- apply(mc$cat,
             c(2),
             quantile,
             probs = c(0.025,0.5,0.975))
cat <- t(cat)
fox <- apply(mc$fox,
             c(2),
             quantile,
             probs = c(0.025,0.5,0.975))
fox <- t(fox)

dat <- rbind(data.frame(group = "Cat",
                        x = data_list$PC1,
                        M = cat[,2],
                        L = cat[,1],
                        H = cat[,3]),
             data.frame(group = "Fox",
                        x = data_list$PC1,
                        M = fox[,2],
                        L = fox[,1],
                        H = fox[,3]))

p5<-ggplot(dat, aes(x=x, y=M, color=group)) +
  #geom_pointrange(aes(ymin=L, ymax=H), show.legend=FALSE) +
  geom_point() +
  geom_smooth(method = 'lm', aes(fill=group), show.legend = FALSE) +
  #scale_color_manual(values = c("#000000")) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="PC1", y="Predator Occupancy", size="", color="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.title=element_blank(), legend.position=c(0.93,0.95), legend.background=element_blank(),
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(ncol=1,byrow=FALSE,
                              keywidth=0.5,
                              keyheight=0.3,
                              default.unit="cm"))
jpeg("./results/suppFig.jpeg", width = 6, height = 4, units = 'in', res = 300)
p5
dev.off()

summary(lm(cat[,2]~data_list$PC1))
summary(lm(fox[,2]~data_list$PC1))
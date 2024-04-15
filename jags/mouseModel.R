model{
  ### PRIORS
  # for initial abundance
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  # for captures
  alpha0 ~ dnorm(0, 0.1)
  alpha1 ~ dnorm(0, 0.1)
  alpha2 ~ dnorm(0, 0.1)
  alpha3 ~ dnorm(0, 0.1) 
  # for persistence (phi)
  phi0 ~ dnorm(0, 0.1)
  phi1 ~ dnorm(0, 0.1)
  phi2 ~ dlogis(0,1)
  phi3 ~ dlogis(0,1)
  phi4 ~ dnorm(0, 0.05)
  phi5 ~ dnorm(0, 0.05)
  phi6 ~ dnorm(0, 0.1)
  # for recruitment (gamma)
  gamma0 ~ dnorm(0, 0.1)
  gamma1 ~ dnorm(0, 0.1)
  # because season is dummy-coded, we have a different coefficient for each season
  for (tmpk in 1:2){
    tmp_gamma[tmpk] ~ dnorm(0, 0.1)
  }
  gamma2[1] <- 0
  for (z in 2:nseason){
    gamma2[z] <- tmp_gamma[z-1]
  }
  # predator probability distributions
  for (j in 1:nsite) {
    logit_cat[j] ~ dnorm(cat.mu[j], cat.tau[j])
    cat.tau[j] <- 1 / cat.sd[j]^2
    cat[j] <- exp(logit_cat[j]) / (1 + exp(logit_cat[j]))
    logit_fox[j] ~ dnorm(fox.mu[j], fox.tau[j])
    fox.tau[j] <- 1 / fox.sd[j]^2
    fox[j] <- exp(logit_fox[j]) / (1 + exp(logit_fox[j]))
  }
  
  ### MODEL
  ## Season 1
  # Latent State (Abundance)
  for (j in 1:nsite) {
    N[j,1] ~ dpois(lambda[j,1])
    log(lambda[j,1]) <- beta0 + beta1*PC1[j]
    # Observation State (Captures)
    for (k in 1:nNight) {
      logit(p[j,k,1]) <- alpha0 + alpha1*moon[j,k,1] + alpha2*jDate[j,k,1] + alpha3*effort[j,k,1]
      y[j,k,1] ~ dbin(p[j,k,1], N[j,1])
    }
    ## Season > 1
    # Latent State
    for (t in 2:nSP) {
      N[j,t] <- P[j,t] + R[j,t]
      P[j,t] ~ dbin(phi[j,t], N[j,t-1])
      logit(phi[j,t]) <- phi0 + phi1*PC1[j] + phi2*cat[j] + phi3*fox[j] + 
        phi4*(cat[j]*PC1[j]) + phi5*(fox[j]*PC1[j]) + phi6*contag[j]
      R[j,t] ~ dpois(gamma[j,t])
      log(gamma[j,t]) <- gamma0 + gamma1*contag[j] + gamma2[season[t]]
      # Observation State
      for (k in 1:nNight) {
        logit(p[j,k,t]) <- alpha0 + alpha1*moon[j,k,t] + alpha2*jDate[j,k,t] + alpha3*effort[j,k,t]
        y[j,k,t] ~ dbin(p[j,k,t], N[j,t])
      }
    }
  }
}

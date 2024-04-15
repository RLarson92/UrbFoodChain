model{
  ## PRIORS
  # latent state priors
  # cat
  a0 ~ dnorm(0,0.25)
  a1 ~ dnorm(0,0.25)
  a2 ~ dnorm(0,0.25)
  # coyote
  b0 ~ dnorm(0,0.25)
  b1 ~ dnorm(0,0.25)
  b2 ~ dnorm(0,0.25)
  # fox
  c0 ~ dnorm(0,0.25)
  c1 ~ dnorm(0,0.25)
  # mink
  d0 ~ dnorm(0,0.25)
  d1 ~ dnorm(0,0.25)
  d2 ~ dnorm(0,0.25)
  # cat/coyote
  e0 ~ dnorm(0,0.25)
  e1 ~ dnorm(0,0.25)
  # cat/fox
  g0 ~ dnorm(0,0.25)
  g1 ~ dnorm(0,0.25)
  # cat/mink
  h0 ~ dnorm(0,0.25)
  h1 ~ dnorm(0,0.25)
  # fox/coyote
  l0 ~ dnorm(0,0.25)
  l1 ~ dnorm(0,0.25)
  # fox/mink
  m0 ~ dnorm(0,0.25)
  m1 ~ dnorm(0,0.25)
  # mink/coyote
  n0 ~ dnorm(0,0.25)
  n1 ~ dnorm(0,0.25)
  
  for(z in 1:nspec){
    p0[z] ~ dnorm(0, 0.25)
    phi[z] ~ dnorm(0, 0.25)
    for (tmpk in 1:3){
      tmp_p[z,tmpk] ~ dnorm(0, 0.25)
    }
    p1[z,1]<-0
    for (k in 2:nseason){
      p1[z,k] <- tmp_p[z,k-1]
    }
  }
  
  ## MODELS
  # loop over sites
  for(j in 1:nsite){
    ### SEASON 1
    ## STATE MODEL - multivariate categorical
    # natural parameters
    # 1 = cat, 2 = coyote, 3 = fox, 4 = mink
    f1[j,1] <- a0 + a1*house[j] + a2*crop[j]
    f2[j,1] <- b0 + b1*prairie[j] + b2*crop[j]
    f3[j,1] <- c0 + c1*forest[j]
    f4[j,1] <- d0 + d1*forest[j] + d2*water[j]
    f12[j,1] <- e0 + e1*imperv[j]
    f13[j,1] <- g0 + g1*imperv[j]
    f14[j,1] <- h0 + h1*imperv[j]
    f23[j,1] <- l0 + l1*imperv[j]
    f24[j,1] <- m0 + m1*imperv[j]
    f34[j,1] <- n0 + n1*imperv[j]
    # Psi gives latent state 'category' for each site
    # for Psi[i,j,k], i = latent occupancy category (e.g., 100, 101, 010), j = individual site, 
    # k = individual season
    # the math below calculates the odds of each latent occupancy category for each iteration
    # all species present
    Psi[1,j,1] <- exp(f1[j,1] + f2[j,1] + f3[j,1] + f12[j,1] + f13[j,1] + f14[j,1] + f23[j,1] + f24[j,1] + f34[j,1])
    # only cats
    Psi[2,j,1] <- exp(f1[j,1])
    # only coyote
    Psi[3,j,1] <- exp(f2[j,1])
    # only fox
    Psi[4,j,1] <- exp(f3[j,1])
    # only mink
    Psi[5,j,1] <- exp(f4[j,1])
    # cat & coyote
    Psi[6,j,1] <- exp(f1[j,1] + f2[j,1] + f12[j,1])
    # cat & fox
    Psi[7,j,1] <- exp(f1[j,1] + f3[j,1] + f13[j,1])
    # cat & mink
    Psi[8,j,1] <- exp(f1[j,1] + f4[j,1] + f14[j,1])
    # coyote & fox
    Psi[9,j,1] <- exp(f2[j,1] + f3[j,1] + f23[j,1])
    # coyote & mink
    Psi[10,j,1] <- exp(f2[j,1] + f4[j,1] + f24[j,1])
    # fox & mink
    Psi[11,j,1] <- exp(f3[j,1] + f4[j,1] + f34[j,1])
    # no species present
    Psi[12,j,1] <- 1
    # model latent occupancy state as categorical random variable derived from categories in Psi
    x[j,1] ~ dcat(Psi[,j,1])
    
    ## OBSERVATION MODEL
    # loop over species
    for(i in 1:nspec){
      # observation modeled as function of trap days (# days camera was active per week per season)
      logit(p[i,j,1]) <- p0[i] + p1[i,season[1]]
      # occupancy modeled as state * observation
      # note: Xcat is a matrix of 2^n rows and n columns, filled w/ 0 & 1, where n = # of species
      # (see above, y is the 0/1/NA array of detections, indexed by species, site, week, season)
        y[i,j,1] ~ dbinom(Xcat[x[j,1],i]*p[i,j,1], J[j,1])
    }
    
    ### SEASON >1
    ## STATE MODEL - multivariate categorical
    # natural parameters
    # for season >1, single species natural parameters have an additional parameter w/ coefficient
    # phi that is multiplied by that site's occupancy state from the previous season
    for(t in 2:nSP){
      f1[j,t] <- a0 + a1*house[j] + a2*crop[j] + phi[1]*Xcat[x[j,t-1],1]
      f2[j,t] <- b0 + b1*prairie[j] + b2*crop[j] + phi[2]*Xcat[x[j,t-1],2]
      f3[j,t] <- c0 + c1*forest[j] + phi[3]*Xcat[x[j,t-1],3]
      f4[j,t] <- d0 + d1*forest[j] + d2*water[j] + phi[4]*Xcat[x[j,t-1],4]
      f12[j,t] <- e0 + e1*imperv[j]
      f13[j,t] <- g0 + g1*imperv[j]
      f14[j,t] <- h0 + h1*imperv[j] 
      f23[j,t] <- l0 + l1*imperv[j]
      f24[j,t] <- m0 + m1*imperv[j]
      f34[j,t] <- n0 + n1*imperv[j]

      Psi[1,j,t] <- exp(f1[j,t] + f2[j,t] + f3[j,t] + f12[j,t] + f13[j,t] + f14[j,t] + f23[j,t] + f24[j,t] + f34[j,t])
      Psi[2,j,t] <- exp(f1[j,t])
      Psi[3,j,t] <- exp(f2[j,t])
      Psi[4,j,t] <- exp(f3[j,t])
      Psi[5,j,t] <- exp(f4[j,t])
      Psi[6,j,t] <- exp(f1[j,t] + f2[j,t] + f12[j,t])
      Psi[7,j,t] <- exp(f1[j,t] + f3[j,t] + f13[j,t])
      Psi[8,j,t] <- exp(f1[j,t] + f4[j,t] + f14[j,t])
      Psi[9,j,t] <- exp(f2[j,t] + f3[j,t] + f23[j,t])
      Psi[10,j,t] <- exp(f2[j,t] + f4[j,t] + f24[j,t])
      Psi[11,j,t] <- exp(f3[j,t] + f4[j,t] + f34[j,t])
      Psi[12,j,t] <- 1

      x[j,t] ~ dcat(Psi[,j,t])
      
      ## OBSERVATION MODEL
      # loop over species
      for(i in 1:nspec){
        logit(p[i,j,t]) <- p0[i] + p1[i,season[t]]
        y[i,j,t] ~ dbinom(Xcat[x[j,t],i]*p[i,j,t], J[j,t])
      }
    }
  }
}


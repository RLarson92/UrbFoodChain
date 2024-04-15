# initial starting values for each chain
rod_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      beta0 = rnorm(1),
      beta1 = rnorm(1),
      alpha0 = rnorm(1),
      alpha1 = rnorm(1),
      alpha2 = rnorm(1),
      alpha3 = rnorm(1),
      phi0 = rnorm(1),
      phi1 = rnorm(1),
      phi2 = rlogis(1),
      phi3 = rlogis(1),
      phi4 = rnorm(1),
      phi5 = rnorm(1),
      phi6 = rnorm(1),
      gamma0 = rnorm(1),
      gamma1 = rnorm(1),
      tmp_gamma = rnorm(2),
      N = N_init,
      R = R_init,
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(
    switch(
      chain,
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
  )
}

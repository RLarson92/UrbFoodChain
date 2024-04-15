# initial starting values for each chain
pred_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      x = matrix(1,
                 nrow = data_list$nsite,
                 ncol = data_list$nSP),
      p0 = rnorm(data_list$nspec), 
      tmp_p = matrix(1,
                     nrow = data_list$nspec,
                     ncol = 3),
      phi = rnorm(data_list$nspec), 
      a0 = rnorm(1),
      a1 = rnorm(1),
      a2 = rnorm(1),
      b0 = rnorm(1),
      b1 = rnorm(1),
      b2 = rnorm(1),
      c0 = rnorm(1),
      c1 = rnorm(1),
      d0 = rnorm(1), 
      d1 = rnorm(1),
      d2 = rnorm(1),
      e0 = rnorm(1), 
      e1 = rnorm(1), 
      g0 = rnorm(1), 
      g1 = rnorm(1),
      h0 = rnorm(1),
      h1 = rnorm(1),
      l0 = rnorm(1),
      l1 = rnorm(1),
      m0 = rnorm(1),
      m1 = rnorm(1),
      n0 = rnorm(1),
      n1 = rnorm(1)
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

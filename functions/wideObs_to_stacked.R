wideObs_to_stacked <- function(input_df, surveys_per_bout){
  nbouts <- ncol(input_df) / surveys_per_bout
  inds <- split(1:(nbouts*surveys_per_bout), rep(1:nbouts, each=surveys_per_bout))
  split_df <- lapply(1:nbouts, function(i){
    out <- input_df[, inds[[i]]]
    out$Site <- moon$Site
    out$Season <- i
    names(out)[1:3] <- paste0("night", 1:3)
    out
  })
  stack_df <- do.call("rbind", split_df)
  stack_df
}
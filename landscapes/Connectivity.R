###########################################################################################################
#                                                                                                         #
#                                 Code for calculating contagion index values                             #
#                                          using `landscapemetrics`                                       #
#                                                                                                         #
#                                         Last edited: 25 Jan 2024                                        #
#                                                                                                         #
###########################################################################################################
# You'll need to load the 'raster' and 'landscapemetrics' packages to perform this analysis
library(raster)
library(landscapemetrics)

# Code to list all the files in the '.tifs' subfolder and import them in the R workspace
rasterList <- list.files("./landscapes/tifs/", full.names = TRUE)
rasterList <- rasterList[-46]  # this is some Windows file not part of the dataset

# telling R that all the files in rasterList are, in fact, rasters
allRasters <- lapply(rasterList, raster)

# creating an empty list item to store the results of the analysis in
resultsList_contag <- list()

# automating the calculation the contagion values for all 45 sites
for(i in 1:length(allRasters)){
  resultsList_contag[[i]] <- lsm_l_contag(allRasters[[i]])
}
# a couple of the trap sites contain only 1 land cover class, and will cause a warning
# from landscapemetrics. When inspected visually in ArcMap or the like, you'll see their
# contagion value should be 0 (i.e., low interspersion and high "clumping"). I hand-edited
# the '0's into my covariate datafile

# writing the results into that blank list and adding the site name
results_contag <- rlist::list.rbind(resultsList_contag)
results_contag$site <- rasterList

# writing the results to a .csv. I was lazy and just copied the values from this .csv
# into the bigger .csv I have that has all the plot covariates, rather than creating
# an R script to stitch them together
write.csv(results_contag, "./data/contag.csv", row.names = FALSE)

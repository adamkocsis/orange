library(tinytest)
library(orange)
library(rgplates) # only for testing setup
library(parallel)

# enforce correct names
if(rgplates:::getOS()=="linux") wd <- file.path(Sys.getenv("Dropbox"), "Software/orange")
if(rgplates:::getOS()=="windows") wd <- file.path("D:/orange")
if(rgplates:::getOS()=="osx") wd <- file.path("~/Desktop/orange")

setwd(wd)

# make a cluster of 8
cl <- parallel::makeCluster(4, outfile="")
parallel::clusterCall(cl, source, "orange/tests/source.R")

# the range metrics
occupancy_results <- run_test_dir("orange/tests/occupancy")
maxdist_results <- run_test_dir("orange/tests/maxdist")
latrange_results <- run_test_dir("orange/tests/latrange")
mstlength_results <- run_test_dir("orange/tests/mstlength")
radius_results <- run_test_dir("orange/tests/radius")

# the shape metrics
centroid_results <- run_test_dir("orange/tests/centroid")
gappiness_results <- run_test_dir("orange/tests/gappiness")

# Small circle tests

# Finish
stopCluster(cl)

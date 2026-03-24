library(tinytest)
library(biodome)
library(rgplates) # only for testing setup
library(parallel)

# enforce correct names
if(rgplates:::getOS()=="linux") wd <- file.path(Sys.getenv("Dropbox"), "Software/biodome")
if(rgplates:::getOS()=="windows") wd <- file.path("D:/biodome")
if(rgplates:::getOS()=="osx") wd <- file.path("~/Desktop/biodome")

setwd(wd)

# make a cluster of 8
cl <- parallel::makeCluster(4, outfile="")
parallel::clusterCall(cl, source, "biodome/tests/source.R")

# the online
occurence <- run_test_dir("biodome/tests/occurrence")

# Finish
stopCluster(cl)

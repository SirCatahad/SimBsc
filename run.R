## Clean up the workspace
rm(list=ls(all=TRUE))

##Libraries
library("SimDesign")

source("init.R")
source("sim_functions.R")

set.seed(103415)

tmp <- simData(n=500, coefficients=coeff_list[[3]], r.sqrd=r.squared[1], parameters=parameters)
tmp <- makeMissing(data=tmp, mechanism="MCAR", percent=.25)
tmp

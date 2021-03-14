## Clean up the workspace
rm(list=ls(all=TRUE))


setwd("C:/Users/wufft/Desktop/RL Kram/Uni/Bachelor Thesis/R files/code/assignment_sim")
##Libraries
library("SimDesign")

source("init.R")
source("sim_functions.R")

set.seed(103415)


start_time <- Sys.time()
dat <- simData(n=500, coefficients=coeff_list[[1]], r.sqrd=r.squared[1], parameters=parameters)
dat <- makeMissing(data=dat, mechanism="MCAR", percent=.25)
end_time <- Sys.time()

sum(is.na(dat))
end_time - start_time

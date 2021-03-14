## Clean up the workspace
rm(list=ls(all=TRUE))

##Libraries
library("SimDesign")

source("init.R")
source("sim_functions.R")

set.seed(103415)


start_time <- Sys.time()
for(i in 1:1000)
{
  dat <- simData(n=500, coefficients=coeff_list[[3]], r.sqrd=r.squared[1], parameters=parameters)
  dat <- makeMissing(data=dat, mechanism="MCAR", percent=.25)
}
end_time <- Sys.time()

sum(is.na(dat))
end_time - start_time

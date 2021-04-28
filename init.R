#Parameters of the simulation:
parameters <- list()



#Variable
covariances <- c(0, 3.33, 10)       # Covariances
r.squared   <- c(0, .1, .5)         # R-square
snr         <- c(.325, .6)          # SNR
mtype       <- c(1,2)               # Matching type for PMM
km          <- c(1, 3 ,5, 10)       # k for each match type
method      <- c("MI", "PMM")

nya <- expand.grid(covariances,snr, km, mtype)

iter <- 10                          # iterations

#Fixed
parameters$n    <- 500              #sample size
parameters$miss <- .25              #percentage of missingness
parameters$M    <- 10               #number of imputations




#Parameters of the simulation:
parameters <- list()

#Variable
betav <- c(0, 3.33, 10)       # Covariances
snr         <- c(.325, .6)          # SNR
mtype       <- c(1,2)               # Matching type for PMM
km          <- c(1, 3 ,5, 10)       # k for each match type
iter <- 1000                        # iterations

#Fixed
parameters$n    <- 500              #sample size
parameters$miss <- .25              #percentage of missingness
parameters$M    <- 10               #number of imputations

study <- "study2"                   # Use "study 2" for quadratic misspecified model
outputDir <- paste0("../data/", study, "/")

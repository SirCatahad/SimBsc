#Parameters of the simulation:
parameters <- list()



#Variable
covariances <- c(0, 3.33, 10)       # Covariances
r.squared   <- c(0, .1, .5)         # R-square
mechanism   <- c("MCAR", "MAR")     # Missing data mechanism
snr         <- c(.65, .75)          # SNR
mtype       <- c(1,2)               # Matching type for PMM
km          <- c(1, 3 ,5, 10)       # k for each match type
method      <- c("MI", "PMM")


iter <- 1                           # iterations

#Fixed
parameters$n    <- 500              #sample size
parameters$miss <- .25              #percentage of missingness
parameters$M    <- 10               #number of imputations
parameters$c    <- list(c(1, 1.5),  #Coefficients for study 1
                        c(0, 2))    #and 2 respectively


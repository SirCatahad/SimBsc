#Parameters of the simulation:
parameters <- list()

#Variable
M <- 10                                     # Number of iterations
nvec <- seq(100,10000,100)                  # sample sizes
r.squared <- c(.1, .33, .5, .66, .9)        # R-squared
coeff_list <- list(c(1, .5, 3.3),
                   c(1, .5, 3.3, 2),
                   c(1, .5, 3.3, 2, 1.5))   # Coefficients for the predictor variables, also number of predictors (length -1)

#Fixed
parameters$covariance <- c(.1)                              #covariances



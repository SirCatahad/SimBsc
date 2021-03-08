# Functions
##-------------------------------------------------------------------------------------------------------------------##
simData <- function(n, r.sqrd, coefficients ,parameters)
{
  #Construct covariance matrix
  nrpred <- length(coefficients) - 1
  sigma <- matrix(parameters$covariance, nrpred, nrpred)
  diag(sigma) <- 1.0
  
  #Generate data
  X <- rmvnorm(n = n, mean = rep(0, nrpred), sigma = sigma)

  #Generate error termn
  beta <- coefficients[-1]
  var.model <- t(beta) %*% cov(X) %*% beta
  var.residual <- (var.model/r.sqrd) - var.model
  U = rnorm(n, mean = 0, sd = sqrt(var.residual))
  
  #compute Y
  Y <- coefficients[1] + X  %*% beta + U 
  data <- data.frame(X)
  
  #Add to data frame
  data$Y <- Y
  
  #Return data
  data
}

##-------------------------------------------------------------------------------------------------------------------##

makeMissing <- function()
{
  #Do stuff
}

##-------------------------------------------------------------------------------------------------------------------##

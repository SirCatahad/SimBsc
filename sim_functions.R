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
# data       - the data frame which should get missing observations
# mechanism  - the mechanism of missing data, by default MCAR
# percent    - the proportion of observations that should be set to missing (NA)
# indexRange - A vector of indices restricting which columns should contain missing values
makeMissing <- function(data, mechanism="MCAR", percent, indices)
{
  df <- data
  if(mechanism=="MAR")
  {
    
  }
  
  #MCAR missing data mechanism
  if(mechanism=="MCAR")
  {
    
    #delete the specified amount of cases in each of the specified columns
    for(i in indices)
    {
      df[sample(1:nrow(data), nrow(df)*percent) , i] <- NA
    }

    #return
    df
  }
  else
  {
    stop("Undefined or unsupported missing data mechanism.")
  }
}

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
# 
#
#
#

# data       - the data frame which should get missing observations
# mechanism  - the mechanism of missing data, by default MCAR
# percent    - the proportion of observations that should be set to missing (NA)
# indexRange - A vector of indices restricting which columns should contain missing values, by default everything except the last.
makeMissing <- function(data, mechanism="MCAR", percent, indexRange=-length(data))
{
  #assign the columns that should contain missing data
  df <- data[,indexRange]
  #size of the new data.frame
  size <- ncol(df) * nrow(df)
  #the amount of observations that need to be deleted in order to obtain the set percentage.
  nrdelete <- size*percent
  
  
  if(mechanism=="MAR")
  {
    
  }
  
  #MCAR missing data mechanism
  #Loop and sample a random observations until the amount of missing
  #values is achieved. If the sampled observation has already been set to zero, try again and don't
  #increase the counter.
  if(mechanism=="MCAR")
  {
    counter <- 0
    while(counter < nrdelete)
    {
      r <- sample(1:nrow(df), 1)
      c <- sample(1:ncol(df), 1)
      paste0("Row: ", r, ", Column: ", c)
      if(!(is.na(df[r,c])))
      {
        df[r,c] <- NA
        counter <- counter + 1
      }
    }
    #replace the columns in the original data set with the ones that now contain missing values
    df <- replaceColumns(data, df)
    df
  }
  else
  {
    #return error
  }
}

##-------------------------------------------------------------------------------------------------------------------##

replaceColumns <- function(original, toreplace)
{
  #Createa vector containing the names of the overlap between the replacement and the original data frame
  change_columns <- colnames(original)[colnames(original) %in% colnames(toreplace)]
  
  #Now replace each column
  for(i in 1:length(change_columns))
  {
    original[change_columns[i]] <- toreplace[change_columns[i]]
  }
  #And return the merged data frame
  return(original)
}

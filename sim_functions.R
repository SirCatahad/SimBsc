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
    #Remember the structure so we can put it back together
    rows <- nrow(df)
    cols <- ncol(df)
    cnames <- colnames(df)
    structure <- matrix(nrow=rows, ncol=cols)

    #Make a vector out of the data frame so we can use sample()
    #replace the values
    some_vector <- unlist(df, use.names = FALSE)

    #creating another vector for the indices is actually faster
    #than unlisting with names
    another_vector <- seq(1,length(some_vector),1)
    tmp <- sample(another_vector, to_delete)



    #Sample from named vector
    #tmp <- sample(some_vector, to_delete)
    some_vector[tmp] <- NA




    #Set the sampled values to NA
    some_vector[names(tmp)] <- NA

    #Restructure the vector to a data frame
    df <- relist(some_vector, structure)
    df <- data.frame(df)

    #Restore original column names
    colnames(df) <- cnames

    #replace the columns in the original data

    df <- replaceColumns(data, df)
    
    #return
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

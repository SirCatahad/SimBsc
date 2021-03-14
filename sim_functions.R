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
# indexRange - A vector of indices restricting which columns should contain missing values, by default everything except the last.
makeMissing <- function(data, mechanism="MCAR", percent, indexRange=-length(data))
{
  #assign the columns that should contain missing data
  #drop=FALSE prevents the drop in dimension if  the number of columns is reduced to 1
  df <- data[,indexRange, drop=FALSE]
  #size of the new data.frame
  size <- ncol(df) * nrow(df)
  #the amount of observations that need to be deleted in order to obtain the set percentage.
  to_delete <- size*percent
  
  
  if(mechanism=="MAR")
  {
    
  }
  
  #MCAR missing data mechanism
  if(mechanism=="MCAR")
  {
    #create a matrix containing row and column indices
    indices <- as.matrix(expand.grid(1:nrow(df), 1:ncol(df)))
    
    #randomly select indices until we have enough
    selected <- indices[sample(nrow(indices), to_delete), ]
    
    #set selected observations to missing
    df[selected] <- NA
    
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

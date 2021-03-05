# Functions
##-------------------------------------------------------------------------------------------------------------------##
simCycle <- function(n, r.sqrd, coefficients ,parameters)
{
  #Construct covariance matrix
  nrpred <- length(coefficients) - 1
  sigma <- matrix(parameters$covariance, nrpred, nrpred)
  diag(sigma) <- 1.0
  
  #Generate data
  X <- rmvnorm(n = n, mean = rep(0,nrpred), sigma = sigma)

  #   t(beta) * cov(X) * beta
  #=> (ß1*cov(X1,X1) + ß2*cov(X1,X2)    ß1*cov(X1,X2)+ß2*cov(X2,X2)) t(ß1 ß2)
  #=> ß1(ß1*cov(X1,X1) + ß2*cov(X1,X2)) + ß2(ß1*cov(X1,X2) + ß2*cov(X2,X2))
  #=> ß1^2 * cov(X1,X1) + ß1*ß2*cov(X1,X2) + ß1*ß2*cov(X1,X2) + ß2^2 *cov(X2,X2)
  #=> ß1^2 * cov(X1,X1) + 2*ß1*ß2*cov(X1,X2) + ß2^2 * cov(X2,X2)
  #=> ß1^2 * Var(X1) + ß2^2*Var(X2) + 2*ß1*ß2*cov(X1,X2)
  #=> Var(ß1*X1 + ß2*X2) => Variance of the regression part.
  #
  # R^2 is the proportion of variance in the dependent variable that is predictable from the independent variables.
  #=> R^2 = Var(model)/var(model)+var(residuals)                    | *var(model)+var(residuals)
  #=> R^2 * var(model) + R^2 * var(residuals) = var(model)          | /R^2
  #=> var(model) + var(residuals) = var(model)/R^2                  | -var(model)
  #=> var(residuals) = (var(model)/R^2) - var(model) 
  
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

runreg <- function(dat)
{
  #run regression model
  reg <- lm(Y ~ ., data = dat)
  
  #build output
  out <- list()
  out$coefficients <- reg$coefficients
  out$residuals <- reg$residuals
  out$r.squared <- summary(reg)$r.squared
  
  out
}

##-------------------------------------------------------------------------------------------------------------------##

drawplots <- function(coefficients, data)
{
  #coefficients <- coeff_list[[1]][-1]
  #data <- l[[1]]
  # take mean slopes for each condition
  mean_slope_DT <- data.frame(matrix(NA, nrow=nrow(data)/M ,ncol=ncol(data)))
  sd_slope_DT <- data.frame(matrix(NA, nrow=nrow(data)/M ,ncol=ncol(data)))
  for(i in 1:(nrow(data)/M))
  {
    a <- ((i-1)*M)+1
    b <- i*M
    
    mean_slope_DT[i,] <- apply(X=data[a:b,],2, FUN=mean)
    sd_slope_DT[i,] <- apply(X=data[a:b,], 2,FUN=sd)
  }
  par(mfrow=c(2,2))
  for(j in 1:ncol(mean_slope_DT))
  {

    plot(mean_slope_DT[,j])
    abline(h=coefficients[j], col="blue")
    plot(sd_slope_DT[,j])
  }
}






# ## Visual inspection
# ##-------------------------------------------------------------------------------------------------------------------##
# par(mfrow=c(3,3))
# #Consistency
# 
# # An OLS estimator is said to be consistent if the estimated values approach the true population parameters as the
# # sample size increases. So we calculate the standard deviation for each sample size:
# sdvec <- data.frame(matrix(nrow=length(nvec), ncol=length(param$mean))-1)
# rownames(sdvec) <- nvec
# for(ssize in 1:length(nvec))
# {
#   #since we know they are ordered by sample size we just have to cut out the correct slice of rows each time
#   sdvec[ssize,] <- apply(slope_DT[ssize:(ssize+M),],2,sd)
# }
# 
# for(c in 1:ncol(sdvec))
# {
#   plot(sdvec[,c], col="red")
#   #abline(h=param$beta[c+1], col="blue")
# }
# #standard error is decreasing with increasing sample size.
# 
# ##-------------------------------------------------------------------------------------------------------------------##
# #Unbiasedness
# 
# 
# #Histogram for the intercepts
# hist(intercept_DT)
# abline(v=param$beta[1], col="blue")
# 
# #Histogram for all slope parameters. Note that the blue line denotes the true population parameter.
# for(nya in 1:length(param$mean))
# {
#   hist(slope_DT[,nya])
#   abline(v=param$beta[nya+1], col="blue")
# }
# # Visually this looks fine! Let's take a look at the numbers:
# # For this, we calculate differences between the true population parameters and the mean of the estimated parameters
# # Also add standard deviation, while we're at it.
# 
# test <- data.frame(NA,matrix(nrow=length(param$mean), ncol=1))
# names(test) <- c("diff" , "SD")
# for(w in 1:length(param$mean))
# {
#   test[w,] <- c(param$beta[w+1] - mean(slope_DT[,w]), sd(slope_DT[,w]))
# }
# ## Let's see!
# test
# 
# #The differences seem very small, as does the standard deviation. This means
# #the estimator is close to the true population parameters and thus unbiased (Probably).
# 
# ##-------------------------------------------------------------------------------------------------------------------##
# 
# #Efficiency
# 
# 
# 
# ##-------------------------------------------------------------------------------------------------------------------##






#######################
# Trash
# Y <- param$beta[1]
# for(i in 1:length(param$mean))
# {
#   Y <- Y + param$beta[i+1]*data[,i]
# }
# Y <- Y + U
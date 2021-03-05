## Clean up the workspace
rm(list=ls(all=TRUE))
setwd("C:/Users/wufft/Desktop/RL Kram/Uni/Bachelor Thesis/R files/code/assignment_sim")
##Libraries
library("SimDesign")
ptm <- proc.time()
source("init.R")
source("sim_functions.R")

#Assignment:
#1) Write a Monte Carlo simulation to demonstrate that the OLS estimates of linear regression slopes are unbiased, 
#   efficient, and consistent.
#2) Extend the simulation you wrote for (1) to show that the results still hold when varying the number of predictor variables 
#   in the model and the size of the R-Squared.

#Storage
slope_DT_rsq <- matrix(0, M*length(nvec)*length(r.squared) , length(coeff_list[[1]])-1)
#intercept_DT_rsq <- rep(0, M*length(nvec)*length(r.squared))


## Generate data for different r-squared
##-------------------------------------------------------------------------------------------------------------------##

for(r in 1:length(r.squared))#
{
  for(n in 1:length(nvec))
  {
    for(i in 1:M)
    {
      #simulate some data
      tmp <- simCycle(n=nvec[n] , coefficients=coeff_list[[1]] , r.sqrd=r.squared[r], parameters=parameters)
      
      #run regression and save intercept/slope
      reg <- runreg(tmp)
      
      # store intercepts and slopes ordered by sample size and condition
      slope_DT_rsq[i+M*(n-1)+M*length(nvec)*(r-1),] <- reg$coefficients[-1]
      #intercept_DT_rsq[i+M*(n-1)+M*length(nvec)*(r-1)] <- reg$coefficients[1]
    }
  }
}

##-------------------------------------------------------------------------------------------------------------------##

#Now analyze the effects of varying r. This function takes the mean of every 10 data points
#So we should see a repeating pattern every (length(nvec)*M)/10 units on the x-axis
drawplots(coefficients=coeff_list[[1]][-1],data=slope_DT_rsq)
#And indeed we do, this shows that the slope estimates approximate the true population parameter (blue line)
#when n increases, even when R^2 changes.

# # An OLS estimator is said to be consistent if the estimated values approach the true population parameters as the
# # sample size increases. As we have seen on the SD plots, this appears to be the case

##-------------------------------------------------------------------------------------------------------------------##

# Biased?
# #Histogram for all slope parameters in each condition. Note that the blue line denotes the true population parameter.

for(k in 1:length(r.squared))
{
  for(nya in 2:length(coeff_list[[1]])-1)
  {
    #Make a histogram for each condition and each slope parameter
    a <- ((k-1)*length(nvec)*M)+1
    b <- k*M*length(nvec)
    hist(slope_DT_rsq[a:b,nya], xlab=paste0("R-squared: ", r.squared[k]), main=paste0("Estimates of slope parameter ", coeff_list[[1]][nya+1]))
    abline(v=coeff_list[[1]][nya+1], col="blue")
  }
}

##-------------------------------------------------------------------------------------------------------------------##

MSE <- vector("list",length(coeff_list[[1]]-1)* length(r.squared) )
#Calculate MSE for efficiency:
for(i in 1:length(r.squared))
{
  for(y in 2:length(coeff_list[[1]])-1)
  {
    a <- i+(i-1)*length(nvec)*M+1
    b <- i*M*length(nvec)
    MSE[[i+(y-1)*length(r.squared)]] <- ((slope_DT_rsq[a:b,y] - (coeff_list[[1]][y]))^2)/length(slope_DT_rsq[a:b,y])
  }
  plot(MSE[[i]])
}
##-------------------------------------------------------------------------------------------------------------------##


## Generate data for different number of predictors


##-------------------------------------------------------------------------------------------------------------------##
l <- vector("list", length(coeff_list))

for(c in 1:length(coeff_list))
{
  slope_DT_coeff <- matrix(0, M*length(nvec), length(coeff_list[[c]])-1)
  #intercept_DT_coeff <- rep(0, M*length(nvec)*length(r.squared))
  for(n in 1:length(nvec))
  {
    for(i in 1:M)
    {
      #simulate some data
      tmp <- simCycle(n=nvec[n] , coefficients=coeff_list[[c]] , r.sqrd=r.squared[1], parameters=parameters)
      
      #run regression and save intercept/slope
      reg <- runreg(tmp)
      
      # store slopes ordered by sample size and condition
      slope_DT_coeff[i+M*(n-1),] <- reg$coefficients[-1]
    }
  }
  l[[c]] <- slope_DT_coeff
}

##-------------------------------------------------------------------------------------------------------------------##

#Now analyze the effects of varying number of predictors
for( a in 1:length(l))
{
  drawplots(coefficients=coeff_list[[a]][-1],data=l[[a]])
}

# # An OLS estimator is said to be consistent if the estimated values approach the true population parameters as the
# # sample size increases. As we have seen on the SD plots, this appears to be the case

##-------------------------------------------------------------------------------------------------------------------##

# Biased?
# #Histogram for all slope parameters in each condition. Note that the blue line denotes the true population parameter.
for(k in 1:length(l))
{
  for(y in 2:length(coeff_list[[k]])-1)
  {
    hist(l[[k]][,y], xlab=paste0("Number of predictors: ", length(coeff_list[[k]])-1), main=paste0("Estimates of slope parameter: ", coeff_list[[k]][y+1]))
    abline(v=coeff_list[[k]][y+1], col="blue")
  }
}

##-------------------------------------------------------------------------------------------------------------------##

# Efficient?

MSE <- vector("list",length(l))
#Calculate MSE for efficiency:
for(i in 1:length(l))
{
  for(y in 2:length(coeff_list[[i]])-1)
  {
    MSE[[i]] <- ((l[[i]][,y] - coeff_list[[k]][y])^2)/length(l[[i]][,y])
  }
  plot(MSE[[i]])
}

##-------------------------------------------------------------------------------------------------------------------##

proc.time() - ptm
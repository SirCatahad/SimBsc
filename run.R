## Clean up the workspace
rm(list=ls(all=TRUE))
setwd("C:/Users/wufft/Desktop/RL Kram/Uni/Bachelor Thesis/R files/code")

## Libraries
library("mice")
library("ggplot2")
library("gridExtra")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(103415)

#Data generation and storage
storage <- data.frame(matrix(0,
                             nrow= 93, #nr of conditions.
                             ncol=  9  #saved variables
                             ))
## Storing all important variables plus an identifier so we can easily order by conditions later on
colnames(storage) <- c("Intercept", 
                       "Slope", 
                       "se_intercept", 
                       "se_slope",
                       "conf.intercept.low",
                       "conf.slope.low",
                       "conf.intercept.high",
                       "conf.slope.high",
                       "Condition")

#Stores the resulting data frame of each iteration so we can combine them later
#Will be a few Megabyte. Probably not a good idea if it gets larger
storage_i <- vector(mode="list",iter)


#Time it! 
start_time <- Sys.time()



## Main Simulation loop -----------------------------------------------------------------------------------------------------

#Study 1
for(i in 1:iter)
{
  con <- 1
  #Create a data table for each iteration containing all relevant information
  for(cr in zip(covariances, r.squared))
  {
    #Generate data with given parameters
    data <- try(with(parameters, simData(n    = n,
                                 coefficients = c[[1]],
                                 covariance   = cr[[1]],
                                 r.sqrd       = cr[[2]])))
    
    
    #Extract relevant information, coefficients, r.squared, etc.
    result         <- analyze(data)
    storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                        result$regsum$coefficients[c(3,4)],
                        result$c_int[,1],
                        result$c_int[,2],
                        con)
    con <- con+1
    
    ###############################################
    ## This is purely for later cosmetics and does not actually do anything
    storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                        result$regsum$coefficients[c(3,4)],
                        result$c_int[,1],
                        result$c_int[,2],
                        con)
    con <- con+1
    storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                        result$regsum$coefficients[c(3,4)],
                        result$c_int[,1],
                        result$c_int[,2],
                        con)
    con <- con+1
    #################################################
    
    
    
    #Poke holes into the dataset
    MCAR_out <- try(with(parameters, makeMissing(data      = data,
                                                 mechanism = "MCAR",
                                                 pm        = miss
                                )))
    data_MCAR <- data
    data_MCAR[MCAR_out$r,1] <- NA
    
    #extract relevant information from the set
    result <- analyze(data_MCAR)
    
    #Store
    storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                        result$regsum$coefficients[c(3,4)],
                        result$c_int[,1],
                        result$c_int[,2],
                        con)
    con <- con+1
    
    dList <- list(data_MCAR)
    for(s in snr)
    {
      MAR_out <- try(with(parameters, makeMissing(data      = data,
                                                  mechanism = "MAR",
                                                  pm        = miss,
                                                  preds     = colnames(data),
                                                  snr       = s
      )))
      data_MAR <- data
      data_MAR[MAR_out$r,1] <- NA
      dList[[length(dList)+1]] <- data_MAR
      
      result <- analyze(data_MAR)
      
      #Store
      storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                          result$regsum$coefficients[c(3,4)],
                          result$c_int[,1],
                          result$c_int[,2],
                          con)
      con <- con+1
      
    }
    
    ## impute with PDMI
    for(d in dList)
    {
      PDMI <- try(mice(data   = d,
                       m      = 10,
                       method = "norm",
                       printFlag = FALSE
      ))
      pdmi_Data <- complete(PDMI, method="long")
      # fit <- with(PDMI, lm(Y ~ X))
      # result <- summary(pool(fit), conf.int=TRUE, conf.level = .9)
      
      #Store
      result <- analyze(pdmi_Data)
      
      #Store
      storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                          result$regsum$coefficients[c(3,4)],
                          result$c_int[,1],
                          result$c_int[,2],
                          con)
      
      # storage[con, ] <- c(result$estimate, 
      #                     result$std.error,
      #                     result$`5 %`,
      #                     result$`95 %`,
      #                     con)
      con <- con+1
    }
    
    ## Impute with PMM
    for(match in mtype)
    {
      for(k in km)
      {
        for(d in dList)
        {
          #PMM <- try(mice.impute.pmm(y         = data_MCAR[,1],
          #                           ry        = MCAR_out$r,
          #                           x         = data_MCAR[,2],
          #                           donors    = k,
          #                           matchtype = match))
          PMM <- try(mice(data   = d,
                          m      = 10,
                          method = "pmm",
                          donors = k,
                          printFlag = FALSE
                          ))
          pmm_Data <- complete(PMM, method="long")
          # fit <- with(PMM, lm(Y ~ X+ I(X^2)))
          # result <- summary(pool(fit), conf.int = TRUE, conf.level = .9)
          # #Store
          # storage[con, ] <- c(result$estimate, 
          #                     result$std.error,
          #                     result$`5 %`,
          #                     result$`95 %`,
          #                     con)
          
          #Store
          result <- analyze(pmm_Data)
          
          #Store
          storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                              result$regsum$coefficients[c(3,4)],
                              result$c_int[,1],
                              result$c_int[,2],
                              con)
          con <- con+1
        }
      }
    }
  }
  #Write out!
  storage_i[[i]] <- storage
  write.table(data.frame(storage), 
        file     = paste0("../data/data_iteration_",i,".csv"))
}
end_time <- Sys.time()

end_time - start_time



df <- do.call("rbind", storage_i)
df <- df[order(df$Condition),]


#About 15 seconds per  iteration
#15 minutes for 60 iterations
#Probably about 4-5 hours for one study
#So 8-10 hours total runtime

## Summarize


## Bias estimations -----------------------------------------------------------------------------------------------------


bias_df <- df
## First, 
bias_df[(bias_df["Condition"]>33),"Intercept"] <- bias_df[(bias_df["Condition"]>33),"Intercept"] - parameters$c[[1]][1]
bias_df[(bias_df["Condition"]>33),"Slope"] <- bias_df[(bias_df["Condition"]>33),"Slope"] - parameters$c[[1]][2]

test <- data.frame(matrix(0,nrow=99,0))

#Predicted value over all iterations for mean value of the intercept
test[c("Condition", "Intercept")] <- aggregate(bias_df["Intercept"], by=bias_df["Condition"], FUN=mean)

#aggregate standard error: SD(\beta_hat)
test[c("Condition" , "Intercept SE")] <- aggregate(bias_df["Intercept"], by=bias_df["Condition"], FUN=sd)

##Calculate confidence intervals
## Confidence interval for \beta_i:
## CI_Beta_i_.95 = [\beta_hat_i - 1.96 x SE(\beta_hat_i), \beta_hat_i + 1.96 x SE(\beta_hat_i)]
test["int.low"] <- sapply(test["Intercept"], FUN= (function(x,y) x - 1.96*y), test["Intercept SE"])
test["int.high"] <- sapply(test["Intercept"], FUN= (function(x,y) x + 1.96*y), test["Intercept SE"])


#Predicted value over all iterations for mean value of the slope
test[c("Condition", "Slope")] <- aggregate(bias_df["Slope"], by=bias_df["Condition"], FUN=mean)

#aggregate standard error: SD(\beta_hat)
test[c("Condition" , "Slope SE")] <- aggregate(bias_df["Slope"], by=bias_df["Condition"], FUN=sd)

#Calculate Monte Carlo confidence intervals
test["slp.low"] <- sapply(test["Slope"], FUN= (function(x,y) x - 1.96*y), test["Slope SE"])
test["slp.high"] <- sapply(test["Slope"], FUN= (function(x,y) x + 1.96*y), test["Slope SE"])




## Coverage of confidence Intervals-----------------------------------------------------------------------------------------------------
## Coverage of confidence intervals: How often does the true population parameter lie in the
## estimated 95% confidence interval?

##New column that tells us whether the true population parameter lies in the confidence interval for each
##Iteration in each condition
df["CI_int_cov"] <- NA

##If the true parameter lies in the confidence interval
df[(df["Condition"]>33), "CI_int_cov"] <- mapply(FUN=(function(x,y) 
                              ifelse((parameters$c[[1]][1] >= x && parameters$c[[1]][1] <= y), TRUE, FALSE)),
                           y=df[(df["Condition"]>33), "conf.intercept.high"],
                           x=df[(df["Condition"]>33), "conf.intercept.low"])

##The true parameter for the \beta = 0 is equal to zero
df[(df["Condition"]<=33), "CI_int_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((0 >= x && 0 <= y), TRUE, FALSE)),
  y=df[(df["Condition"]<=33), "conf.intercept.high"],
  x=df[(df["Condition"]<=33), "conf.intercept.low"])

df["CI_slp_cov"] <- NA

##If the true parameter lies in the confidence interval
df[(df["Condition"]>33), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((parameters$c[[1]][2] >= x && parameters$c[[1]][2] <= y), TRUE, FALSE)),
  y=df[(df["Condition"]>33), "conf.slope.high"],
  x=df[(df["Condition"]>33), "conf.slope.low"])

##The true parameter for the \beta = 0 is equal to zero
df[(df["Condition"]<=33), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((0 >= x && 0 <= y), TRUE, FALSE)),
  y=df[(df["Condition"]<=33), "conf.slope.high"],
  x=df[(df["Condition"]<=33), "conf.slope.low"])


##Create new data.frame
CI_df <- data.frame(matrix(0,nrow=99,0))

##Percentage of CIs that contain the true population parameter
CI_df[c("Condition", "CI_int_cov")] <- aggregate(df$CI_int_cov, by=df["Condition"], FUN= function(x) (sum(x)/iter)*100)
CI_df[c("Condition", "CI_slp_cov")] <- aggregate(df$CI_slp_cov, by=df["Condition"], FUN= function(x) (sum(x)/iter)*100)

## Compute Monte Carlo CIs
CI_df["CIcov.int.low"] <- as.vector(sapply(CI_df["CI_int_cov"], FUN= (function(x) x - 1.96*x/parameters$n)))
CI_df["CIcov.int.high"] <- as.vector(sapply(CI_df["CI_int_cov"], FUN= (function(x) x + 1.96*x/parameters$n)))

## Compute Monte Carlo CIs
CI_df["CIcov.slp.low"] <- as.vector(sapply(CI_df["CI_slp_cov"], FUN= (function(x) x - 1.96*x/parameters$n)))
CI_df["CIcov.slp.high"] <- as.vector(sapply(CI_df["CI_slp_cov"], FUN= (function(x) x + 1.96*x/parameters$n)))


## Empirical standard Error-----------------------------------------------------------------------------------------------------
## The standard deviation of β over 1,000 replications (henceforth the ‘empirical standard error’)

## Empirical SE is calculated as: sqrt(1/(nsim-1) * sum_i=1_nsim ( Theta_hat_i - Theta_mean)^2)

SE_df <- data.frame(matrix(0,nrow=99,0))

df["Slope Mean"] <- ave(x=unlist(df["Slope"]), df["Condition"])


df["Empirical SE Slope"] <- mapply(FUN= function(th, tm) sqrt(sum((th - tm)^2)/(iter-1)), 
                                      df$Slope, df$`Slope Mean`)


sqrt((1/(iter-1)) * sum((unlist(df$Slope) - df$`Slope Mean`)^2))

SE_df <- aggregate(df$Slope, by=df["Condition"], FUN= function(th, tm) sqrt((1/(iter-1)) * sum((th - tm)^2)), df$`Slope Mean`)


## Data prep-----------------------------------------------------------------------------------------------------
## Data method used
methods <- c(rep("Complete",3),
             rep("Incomplete",3),
             rep("PDMI", 3),
             rep("PMM - Type I", 12),
             rep("PMM - Type II", 12))


## Identifier for k's used
kv <- c(rep(" ", 9),
        rep("k =  1",3),
        rep("k =  3",3),
        rep("k =  5",3),
        rep("k = 10",3),
        rep("k =  1",3),
        rep("k =  3",3),
        rep("k =  5",3),
        rep("k = 10",3)
)

test["k"]      <- kv
test$k <- factor(test$k)

test["method"] <- methods

test <- data.frame(test)


vcol <- c(rep("Complete",3),
          rep(c("MCAR","Weak MAR","Strong MAR"),10))

test["color"] <- vcol


combined_df <- merge(test, CI_df, by="Condition")


##Split by covariance for better overview
plot_df <- split(combined_df, rep(1:3, length.out = nrow(test), each = ceiling(nrow(test)/3)))

plotlist <- list()


## Plots, drawings and everything neat -----------------------------------------------------------------------------------------------------
## Bias Intercept Plotgroup
plot1 <- makePlot(data   = plot_df[[1]], 
                  xint   = 0, 
                  title  = "β = 0", 
                  x      = plot_df[[1]]$Intercept, 
                  xlow   = plot_df[[1]]$int.low, 
                  xhigh  = plot_df[[1]]$int.high,
                  y      = plot_df[[1]]$Condition,
                  xlimits= c(-1,1))

plot2 <- makePlot(data   = plot_df[[2]], 
                  xint   = 0, 
                  title  = "β = 3.33", 
                  x      = plot_df[[2]]$Intercept, 
                  xlow   = plot_df[[2]]$int.low, 
                  xhigh  = plot_df[[2]]$int.high,
                  y      = plot_df[[2]]$Condition,
                  xlimits= c(-1.5,1))

plot3 <- makePlot(data   = plot_df[[3]], 
                  xint   = 0, 
                  title  = "β = 10", 
                  x      = plot_df[[3]]$Intercept, 
                  xlow   = plot_df[[3]]$int.low, 
                  xhigh  = plot_df[[3]]$int.high,
                  y      = plot_df[[3]]$Condition,
                  xlimits= c(-1,1))


grid.arrange(plot1, plot2, plot3, ncol=3, nrow=1, top="Bias Intercept")


##Bias Slope plot group
plot2_1 <- makePlot(data   = plot_df[[1]], 
                    xint   = 0, 
                    title  = "β = 0", 
                    x      = plot_df[[1]]$Slope, 
                    xlow   = plot_df[[1]]$slp.low, 
                    xhigh  = plot_df[[1]]$slp.high,
                    y      = plot_df[[1]]$Condition,
                    xlimits= c(-1,1))

plot2_2 <- makePlot(data   = plot_df[[2]], 
                    xint   = 0, 
                    title  = "β = 3.33", 
                    x      = plot_df[[2]]$Slope, 
                    xlow   = plot_df[[2]]$slp.low, 
                    xhigh  = plot_df[[2]]$slp.high,
                    y      = plot_df[[2]]$Condition,
                    xlimits= c(-1.5,1.5))

plot2_3 <- makePlot(data   = plot_df[[3]], 
                    xint   = 0, 
                    title  = "β = 10", 
                    x      = plot_df[[3]]$Slope, 
                    xlow   = plot_df[[3]]$slp.low, 
                    xhigh  = plot_df[[3]]$slp.high,
                    y      = plot_df[[3]]$Condition,
                    xlimits= c(-1,1))


grid.arrange(plot2_1, plot2_2, plot2_3, ncol=3, nrow=1, top="Bias Slope")


## Confidence Interval coverage plot group
plot3_1 <- makePlot(data   = plot_df[[1]], 
                    xint   = 0, 
                    title  = "β = 0", 
                    x      = plot_df[[1]]$CI_int_cov, 
                    xlow   = plot_df[[1]]$CIcov.int.low, 
                    xhigh  = plot_df[[1]]$CIcov.int.high,
                    y      = plot_df[[1]]$Condition,
                    xlimits= c(50,105))

plot3_2 <- makePlot(data   = plot_df[[2]], 
                    xint   = 0, 
                    title  = "β = 3.33", 
                    x      = plot_df[[2]]$CI_int_cov, 
                    xlow   = plot_df[[2]]$CIcov.int.low, 
                    xhigh  = plot_df[[2]]$CIcov.int.high,
                    y      = plot_df[[2]]$Condition,
                    xlimits= c(50,105))

plot3_3 <- makePlot(data   = plot_df[[3]], 
                    xint   = 0, 
                    title  = "β = 10", 
                    x      = plot_df[[3]]$CI_int_cov, 
                    xlow   = plot_df[[3]]$CIcov.int.low, 
                    xhigh  = plot_df[[3]]$CIcov.int.high,
                    y      = plot_df[[3]]$Condition,
                    xlimits= c(50,105))


grid.arrange(plot3_1, plot3_2, plot3_3, ncol=3, nrow=1, top="CI Coverage Intercept")

## Confidence Interval coverage plot group
plot4_1 <- makePlot(data   = plot_df[[1]], 
                    xint   = 0, 
                    title  = "β = 0", 
                    x      = plot_df[[1]]$CI_slp_cov, 
                    xlow   = plot_df[[1]]$CIcov.slp.low, 
                    xhigh  = plot_df[[1]]$CIcov.slp.high,
                    y      = plot_df[[1]]$Condition,
                    xlimits= c(50,105))

plot4_2 <- makePlot(data   = plot_df[[2]], 
                    xint   = 0, 
                    title  = "β = 3.33", 
                    x      = plot_df[[2]]$CI_slp_cov, 
                    xlow   = plot_df[[2]]$CIcov.slp.low, 
                    xhigh  = plot_df[[2]]$CIcov.slp.high,
                    y      = plot_df[[2]]$Condition,
                    xlimits= c(50,105))

plot4_3 <- makePlot(data   = plot_df[[3]], 
                    xint   = 0, 
                    title  = "β = 10", 
                    x      = plot_df[[3]]$CI_slp_cov, 
                    xlow   = plot_df[[3]]$CIcov.slp.low, 
                    xhigh  = plot_df[[3]]$CIcov.slp.high,
                    y      = plot_df[[3]]$Condition,
                    xlimits= c(50,105))


grid.arrange(plot3_1, plot3_2, plot3_3, ncol=3, nrow=1, top="CI Coverage Intercept")

## Empirical standard error plot group



grid.arrange(plot3_1, plot3_2, plot3_3, ncol=3, nrow=1, top="CI Coverage Intercept")




## TO - DO List -----------------------------------------------------------------------------------------------------


## TO - DO
## Find out how exactly MCCIs work                            (  -  )
## Empirical Standard Errors                                  (https://cran.r-project.org/web/packages/rsimsum/vignettes/A-introduction.html)
## cbind all data frames, plot stuff                          (almost, as soon as MCCIs are done)
## Figure out how to properly space in ggplot                 (done)
## fix xlimits                                                (as soon as MCCIs are done)


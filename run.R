## Clean up the workspace
rm(list=ls(all=TRUE))

## Libraries
library("mice")
library("ggplot2")
library("gridExtra")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(152351)

#Data generation and storage
storage <- data.frame(matrix(0,
                             nrow= 99, #nr of conditions.
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
      # result <- summary(pool(fit), conf.int=TRUE)
      
      #Store
      result <- analyze(pdmi_Data)
      
      storage[con, ] <- c(result$regsum$coefficients[c(1,2)], 
                          result$regsum$coefficients[c(3,4)],
                          result$c_int[,1],
                          result$c_int[,2],
                          con)
      
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
## First, zero all conditions
complete_df <- data.frame(matrix(0,nrow=99,0))

## Bias is calculated as: 1/iter * sum_i=1_nsim ( theta_hat_i - theta ).
## First, we get each deviation
bias_df[(bias_df["Condition"]>33),"Slope"] <- bias_df[(bias_df["Condition"]>33),"Slope"] - covariances[2]
bias_df[(bias_df["Condition"]>66),"Slope"] <- bias_df[(bias_df["Condition"]>66),"Slope"] - (covariances[3]-covariances[2])

## Then we sum (per condition) and divide by nsum(iter)
complete_df[c("Condition", "Slope")] <- aggregate(bias_df["Slope"], by=bias_df["Condition"], FUN=(function(x) sum(x)/iter))

## MCSE(Bias) is calculated by: sqrt(1/nsim-1 * sum_i=1_nsim((theta_hat_i - theta_mean)^2)/nsim). This means we need a
## Theta_mean by condition
bias_df["Slope Mean Group"] <- ave(x=unlist(df["Slope"]), df["Condition"])


## Now calculate Monte carlo standard error of bias (MCSE(Bias))
complete_df[c("Condition", "Slope.SE")] <- aggregate(bias_df["Slope"], by=bias_df["Condition"],
                                                     FUN=function(x) sqrt((1/(iter-1)*sum((x^2)))/iter ))

## And finally Monte Carlo CIs
complete_df["slp.low"] <- sapply(complete_df["Slope"], FUN= (function(x,y) x - 1.96*y), complete_df["Slope.SE"])
complete_df["slp.high"] <- sapply(complete_df["Slope"], FUN= (function(x,y) x + 1.96*y), complete_df["Slope.SE"])

rm(bias_df)


## Coverage of confidence Intervals-----------------------------------------------------------------------------------------------------
## Coverage of confidence intervals: How often does the true population parameter lie in the
## estimated 95% confidence interval?

##New column that tells us whether the true population parameter lies in the confidence interval for each
##Iteration in each condition
df["CI_slp_cov"] <- NA

##Logical column whether the true parameter lies in the 95% confidence interval

df[(df["Condition"]>33), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((covariances[2] >= x && covariances[2] <= y), 1, 0)),
  y=df[(df["Condition"]>33), "conf.slope.high"],
  x=df[(df["Condition"]>33), "conf.slope.low"])

df[(df["Condition"]>66), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((covariances[3] >= x && covariances[3] <= y), 1, 0)),
  y=df[(df["Condition"]>66), "conf.slope.high"],
  x=df[(df["Condition"]>66), "conf.slope.low"])

df[(df["Condition"]<=33), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((covariances[1] >= x && covariances[1] <= y), 1, 0)),
  y=df[(df["Condition"]<=33), "conf.slope.high"],
  x=df[(df["Condition"]<=33), "conf.slope.low"])




##Percentage of CIs that contain the true population parameter
complete_df[c("Condition", "CI_slp_cov")] <- aggregate(df$CI_slp_cov, by=df["Condition"], FUN= function(x) (sum(x)))
complete_df["CI_slp_cov"] <- complete_df["CI_slp_cov"]/iter


complete_df["CI_slp_cov"] <- complete_df["CI_slp_cov"]*100


## SE of the slope CI coverage
complete_df[c("Condition", "CI.slp.se")] <- as.vector(sapply(complete_df$CI_slp_cov, FUN=function(x) sqrt((x*(100-x))/iter)))

## Compute Monte Carlo CIs
complete_df["CIcov.slp.low"] <- as.vector(sapply(complete_df["CI_slp_cov"], FUN= (function(x,y) x - 1.96*y), complete_df$CI.slp.se ))
complete_df["CIcov.slp.high"] <- as.vector(sapply(complete_df["CI_slp_cov"], FUN= (function(x,y) x + 1.96*y), complete_df$CI.slp.se ))


## Empirical standard Error-----------------------------------------------------------------------------------------------------
## The standard deviation of β over 1,000 replications (henceforth the ‘empirical standard error’)

## Empirical SE is calculated as: sqrt(1/(nsim-1) * sum_i=1_nsim ( Theta_hat_i - Theta_mean)^2)
df["Slope Mean"] <- ave(x=unlist(df["Slope"]), df["Condition"])


df["Empirical SE Slope"] <- mapply(FUN= function(th, tm) sqrt(sum((th - tm)^2)/(iter-1)), 
                                   df$Slope, 
                                   df$`Slope Mean`)



df["tmp"] <- mapply(FUN=function(th, tm) (th-tm)^2,
                    df$Slope,
                    df$`Slope Mean`)


complete_df[c("Condition","Empirical SE")] <- aggregate(df$tmp, by=df["Condition"], FUN= function(x) sum(x))

complete_df["Empirical SE"] <- as.vector(sapply(complete_df["Empirical SE"], FUN= function(x) sqrt(1/(iter-1) * x )))

#sqrt((1/(iter-1)) * sum((th - tm)^2)), df$`Slope Mean`)

#Standard Error of the standard error. Meta Standard Error. Have I gone too far?
complete_df[c("Condition","SE SE")] <- aggregate(df$`Empirical SE Slope`, by=df["Condition"], FUN=sd)

complete_df["SE.slp.low"] <- as.vector(sapply(complete_df["Empirical SE"], FUN= (function(x,y) x - 1.96*y/iter), complete_df$`SE SE` ))
complete_df["SE.slp.high"] <- as.vector(sapply(complete_df["Empirical SE"], FUN= (function(x,y) x + 1.96*y/iter), complete_df$`SE SE` ))

## Data prep-----------------------------------------------------------------------------------------------------

## Method used per conditions
complete_df["method"] <- c(rep("Complete",3),
                    rep("Incomplete",3),
                    rep("PDMI", 3),
                    rep("PMM - Type I", 12),
                    rep("PMM - Type II", 12))

## Identifier for k's used
complete_df["k"]      <- c(rep(" ", 9),
                    rep("k =  1",3),
                    rep("k =  3",3),
                    rep("k =  5",3),
                    rep("k = 10",3),
                    rep("k =  1",3),
                    rep("k =  3",3),
                    rep("k =  5",3),
                    rep("k = 10",3))

## Add color marker for Conditions
complete_df["color"] <- c(rep("Complete",3),
                   rep(c("MCAR","Weak MAR","Strong MAR"),10))

complete_df <- data.frame(complete_df)


##Split by covariance for better overview
plot_df <- split(complete_df, rep(1:3, length.out = nrow(complete_df), each = ceiling(nrow(complete_df)/3)))

plotlist <- list()


## Plots, drawings and everything neat -----------------------------------------------------------------------------------------------------

##Bias Slope plot group
plot1_1 <- makePlot(data   = plot_df[[1]], 
                    xint   = 0, 
                    title  = "β = 0", 
                    x      = plot_df[[1]]$Slope, 
                    xlow   = plot_df[[1]]$slp.low, 
                    xhigh  = plot_df[[1]]$slp.high,
                    y      = plot_df[[1]]$Condition,
                    xlimits= c(-1,1))

plot1_2 <- makePlot(data   = plot_df[[2]], 
                    xint   = 0, 
                    title  = "β = 3.33", 
                    x      = plot_df[[2]]$Slope, 
                    xlow   = plot_df[[2]]$slp.low, 
                    xhigh  = plot_df[[2]]$slp.high,
                    y      = plot_df[[2]]$Condition,
                    xlimits= c(-1.5,1.5))

plot1_3 <- makePlot(data   = plot_df[[3]], 
                    xint   = 0, 
                    title  = "β = 10", 
                    x      = plot_df[[3]]$Slope, 
                    xlow   = plot_df[[3]]$slp.low, 
                    xhigh  = plot_df[[3]]$slp.high,
                    y      = plot_df[[3]]$Condition,
                    xlimits= c(-1,1))


grid.arrange(plot1_1, plot1_2, plot1_3, ncol=3, nrow=1, top="Bias")

## Confidence Interval coverage plot group
plot2_1 <- makePlot(data   = plot_df[[1]], 
                    xint   = 95, 
                    title  = "β = 0", 
                    x      = plot_df[[1]]$CI_slp_cov, 
                    xlow   = plot_df[[1]]$CIcov.slp.low, 
                    xhigh  = plot_df[[1]]$CIcov.slp.high,
                    y      = plot_df[[1]]$Condition,
                    xlimits= c(50,105))

plot2_2 <- makePlot(data   = plot_df[[2]], 
                    xint   = 95, 
                    title  = "β = 3.33", 
                    x      = plot_df[[2]]$CI_slp_cov, 
                    xlow   = plot_df[[2]]$CIcov.slp.low, 
                    xhigh  = plot_df[[2]]$CIcov.slp.high,
                    y      = plot_df[[2]]$Condition,
                    xlimits= c(50,105))

plot2_3 <- makePlot(data   = plot_df[[3]], 
                    xint   = 95, 
                    title  = "β = 10", 
                    x      = plot_df[[3]]$CI_slp_cov, 
                    xlow   = plot_df[[3]]$CIcov.slp.low, 
                    xhigh  = plot_df[[3]]$CIcov.slp.high,
                    y      = plot_df[[3]]$Condition,
                    xlimits= c(50,105))


grid.arrange(plot2_1, plot2_2, plot2_3, ncol=3, nrow=1, top="Coverage of 95% Confidence Intervals")

## Empirical standard error plot group
plot3_1 <- makePlot(data   = plot_df[[1]], 
                    xint   = 0, 
                    title  = "β = 0", 
                    x      = plot_df[[1]]$Empirical.SE, 
                    xlow   = plot_df[[1]]$SE.slp.low, 
                    xhigh  = plot_df[[1]]$SE.slp.high,
                    y      = plot_df[[1]]$Condition,
                    xlimits= c(0,1))

plot3_2 <- makePlot(data   = plot_df[[2]], 
                    xint   = 0, 
                    title  = "β = 3.33", 
                    x      = plot_df[[2]]$Empirical.SE, 
                    xlow   = plot_df[[2]]$SE.slp.low, 
                    xhigh  = plot_df[[2]]$SE.slp.high,
                    y      = plot_df[[2]]$Condition,
                    xlimits= c(0,1))

plot3_3 <- makePlot(data   = plot_df[[3]], 
                    xint   = 0, 
                    title  = "β = 10", 
                    x      = plot_df[[3]]$Empirical.SE, 
                    xlow   = plot_df[[3]]$SE.slp.low, 
                    xhigh  = plot_df[[3]]$SE.slp.high,
                    y      = plot_df[[3]]$Condition,
                    xlimits= c(0,2))


grid.arrange(plot3_1, plot3_2, plot3_3, ncol=3, nrow=1, top="Empirical Standard Error")





## TO - DO List -----------------------------------------------------------------------------------------------------


## TO - DO
## Take another look at PDMI, something seems weird           (nope)
## Find out how exactly MCCIs work                            (https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086 ?)
## Empirical Standard Errors                                  (https://cran.r-project.org/web/packages/rsimsum/vignettes/A-introduction.html)
## cbind all data frames, plot stuff                          (done)
## Figure out how to properly space in ggplot                 (done)
## fix xlimits                                                (kinda done)


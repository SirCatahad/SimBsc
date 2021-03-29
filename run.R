## Clean up the workspace
rm(list=ls(all=TRUE))

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
          # fit <- with(PMM, lm(Y ~ X))
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

# Subtract true parameters from estimates to get Bias
df[(df["Condition"]>31),"Intercept"] <- df[(df["Condition"]>31),"Intercept"] - parameters$c[[1]][1]
df[(df["Condition"]>31),"Slope"] <- df[(df["Condition"]>31),"Slope"] - parameters$c[[1]][2]



#About 15 seconds per  iteration
#15 minutes for 60 iterations
#Probably about 4-5 hours for one study
#So 8-10 hours total runtime

## Summarize
test <- data.frame(matrix(0,nrow=93,0))

test[c("Condition", "Intercept")] <- aggregate(df["Intercept"], by=df["Condition"], FUN=mean)
tmp <- data.frame(aggregate(df["Intercept"], by=df["Condition"], FUN=quantile, probs=c(.025, .975)))
test[c("int.low", "int.high")] <- c(tmp[2]$Intercept[,1], tmp[2]$Intercept[,2])

test[c("Condition", "Slope")] <- aggregate(df["Slope"], by=df["Condition"], FUN=mean)
tmp <- data.frame(aggregate(df["Slope"], by=df["Condition"], FUN=quantile, probs=c(.025, .975)))
test[c("slp.low", "slp.high")] <- c(tmp[2]$Slope[,1], tmp[2]$Slope[,2])

## Data method used
methods <- c("Complete",
             rep("Incomplete",3),
             rep("PDMI", 3),
             rep("PMM - Type I", 12),
             rep("PMM - Type II", 12))


## Identifier for k's used
kv <- c(rep(" ", 7),
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
test["method"] <- methods


##Split by covariance for better overview
plot_df <- split(test, rep(1:3, length.out = nrow(test), each = ceiling(nrow(test)/3)))

plotlist <- list()


## Add color column (also identifies missing data condition)
vcol <- c("Complete",
          rep(c("MCAR","Weak MAR","Strong MAR"),10))


plot_df$`1`["color"] <- vcol
plot_df$`2`["color"] <- vcol
plot_df$`3`["color"] <- vcol
  




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

plot2_1 <- makePlot(data   = plot_df[[1]], 
                  xint   = 0, 
                  title  = "β = 0", 
                  x      = plot_df[[1]]$Slope, 
                  xlow   = plot_df[[1]]$slp.low, 
                  xhigh  = plot_df[[1]]$slp.high,
                  y      = plot_df[[1]]$Condition,
                  xlimits= c(-.5,.5))

plot2_2 <- makePlot(data   = plot_df[[2]], 
                  xint   = 0, 
                  title  = "β = 3.33", 
                  x      = plot_df[[2]]$Slope, 
                  xlow   = plot_df[[2]]$slp.low, 
                  xhigh  = plot_df[[2]]$slp.high,
                  y      = plot_df[[2]]$Condition,
                  xlimits= c(-1,1))

plot2_3 <- makePlot(data   = plot_df[[3]], 
                  xint   = 0, 
                  title  = "β = 10", 
                  x      = plot_df[[3]]$Slope, 
                  xlow   = plot_df[[3]]$slp.low, 
                  xhigh  = plot_df[[3]]$slp.high,
                  y      = plot_df[[3]]$Condition,
                  xlimits= c(-.3,.3))


grid.arrange(plot2_1, plot2_2, plot2_3, ncol=3, nrow=1, top="Bias Slope")


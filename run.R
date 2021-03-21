## Clean up the workspace
rm(list=ls(all=TRUE))
setwd("C:/Users/wufft/Desktop/RL Kram/Uni/Bachelor Thesis/R files/code")

## Libraries
library("mice")
library("ggplot2")

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
    
    nya <- lm(Y ~ X, data)
    nya <- confint(nya, level=.9)
    
    summary(nya)
    
    
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
    for(d in dList)
    {
      ## Impute with PMM
      for(match in mtype)
      {
        for(k in km) #2*4*3*3 + 4
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
          fit <- with(PMM, lm(Y ~ X))
          result <- summary(pool(fit), conf.int = TRUE, conf.level = .9)
          #Store
          storage[con, ] <- c(result$estimate, 
                              result$std.error,
                              result$`5 %`,
                              result$`95 %`,
                              con)
          con <- con+1
        }
      }
      ## impute with PDMI
      PDMI <- try(mice(data   = d,
                      m      = 10,
                      method = "norm",
                      printFlag = FALSE
                      ))
      pdmi_Data <- complete(PDMI, method="long")
      fit <- with(PDMI, lm(Y ~ X))
      result <- summary(pool(fit), conf.int=TRUE, conf.level = .9)
      
      #Store
      storage[con, ] <- c(result$estimate, 
                          result$std.error,
                          result$`5 %`,
                          result$`95 %`,
                          con)
      con <- con+1
    }
  }
  #Write out!
  storage_i[[i]] <- storage
  write.table(data.frame(storage), 
        file     = paste0("../data/data_iteration_",i,".csv"))
}

df <- do.call("rbind", storage_i)
df <- super_df[order(super_df$Condition),]

#Bias, coverage of confidence intervals, and a measure of (in-)efficiency, 
#the standard deviation of β over 1,000 replications (henceforth the ‘empirical standard error’), are summarised
# - Bias
# - coverage of confidence intervals
# - measure of (in-)efficiency
# - standard error

ggplot(super_df)



end_time <- Sys.time()

end_time - start_time
#About 15 seconds per  iteration
#15 minutes for 60 iterations
#Probably about 4-5 hours for one study
#So 8-10 hours total runtime


test <- data.frame(aggregate(df, by=df["Condition"], FUN=mean))



ggplot(test, aes(x=-Condition, ymin=conf.intercept.low, y=Intercept, ymax=conf.intercept.high))+
  geom_pointrange()+
  #geom_linerange(aes(x=Intercept, ymin=Condition, ymax=Condition, xmin=conf.intercept.low, xmax=conf.intercept.high)) +
  labs(title="Coefficients plot")+
  coord_flip()

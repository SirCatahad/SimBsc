## Clean up the workspace
rm(list=ls(all=TRUE))

## Libraries
library("mice")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(162551)

#Stores the resulting data frame of each iteration so we can combine them later
#Will be a few Megabytes. Probably not a good idea if it gets larger
storage_i <- vector(mode="list",iter)

#Time it! 
start_time <- Sys.time()

## Main Simulation loop -----------------------------------------------------------------------------------------------------
#Study 1
for(i in 1:iter)
{
  ## Run one iteration
  storage_i[[i]] <- runIteration(covariances, parameters, km, mtype, snr, "study2", i)
  
  ##Save to disk every 50 iterations
  if(i %% 50==0)
  {
    saveRDS(storage_i[!sapply(storage_i,is.null)],
            file = paste0("../data/study2/data_iteration_",i-49,"_",i,".rds")
    )
    storage_i <- vector(mode="list",iter)
  }
  #  savetodisk
}
end_time <- Sys.time()

end_time - start_time


## Read in all the files
files <- list.files(path = "../data/study2/", pattern = "\\.rds$", full.names = TRUE)
df <- lapply(files, readRDS)

## We now have a list of lists of data frames
df <- do.call("rbind", df)
##After this we have a list of data frames
df <- do.call("rbind", df)
##Now we have combined data frame

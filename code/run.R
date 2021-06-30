## Clean up the workspace
rm(list=ls(all=TRUE))

## Libraries
library("mice")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(162551)

#Stores the resulting data frame of each iteration so we can combine them later
storage_i <- vector(mode="list",50)

## Main Simulation loop -----------------------------------------------------------------------------------------------------
for(i in 1:iter)
{
  ## Run one iteration and save in list
  storage_i[[((i-1)%%50)+1]] <- runIteration(betav, parameters, km, mtype, snr, study, i)
  
  ##Save to disk every 50 iterations
  if(i %% 50==0)
  {
    saveRDS(storage_i[!sapply(storage_i,is.null)],
            file = paste0(outputDir,"data_iteration_",i-49,"_",i,".rds") 
    )
    storage_i <- vector(mode="list",50)
  }
}

## Read in all the files
files <- list.files(path = outputDir, pattern = "\\.rds$", full.names = TRUE)
df <- lapply(files, readRDS)

## We now have a list of lists of data frames
df <- do.call("rbind", df)
##After this we have a list of data frames
df <- do.call("rbind", df)
##Now we have the combined data frame
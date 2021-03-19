## Clean up the workspace
rm(list=ls(all=TRUE))


## Libraries
library("mice")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(103415)

#Data generation and storage

for(i in iter)
{
  #Create a data table for each iteration containing all relevant information
  table <- list()
  for(cr in zip(covariances, r.squared))
  {
    #Generate data with given parameters
    data <- try(with(parameters, simData(n    = n,
                                 coefficients = c[[1]],
                                 covariance   = cr[[1]],
                                 r.sqrd       = cr[[2]])))
    
    
    #Extract relevant information, coefficients, r.squared, etc.
    #analyze(data)
    
    #Poke holes into the dataset
    for(m in mechanism)
    {
      #If it's MCAR
      if(m == "MCAR")
      {
        MCAR_out <- try(with(parameters, makeMissing(data      = data,
                                                     mechanism = m,
                                                     pm        = miss
                                    )))
        data_MCAR <- data
        data_MCAR[MCAR_out$r,1] <- NA
        
        #extract relevant information from the set
        #analyze(data_MCAR)
        
        # Impute
        for(me in method)
        {
          if(me=="PMM")
          {
            for(match in mtype)
            {
              for(k in km)
              {
                
              }
            }
          }
          else
          {
            
          }
        }
      }
      #Otherwise it's MAR
      else
      {
        data_MAR <- data
        for(s in snr)
        {
          MAR_out <- try(with(parameters, makeMissing(data      = data,
                                                      mechanism = m,
                                                      pm        = miss,
                                                      preds     = colnames(data),
                                                      snr       = s
                                                      )))
          data_MAR <- data
          data_MAR[MAR_out$r,1] <- NA
          
          #Impute
          for(me in method)
          {
            if(me=="PMM")
            {
              for(match in mtype)
              {
                for(k in km)
                {
                  #Analyze and store relevant information
                }
              }
            }
            else
            {
              #analyze and store relevant information
            }
          }
        }
      }
    }
  }
}

plot(data)
points(data_MCAR, col='white')
points(data_MAR, col='white')

## Clean up the workspace
rm(list=ls(all=TRUE))


## Libraries
library("mice")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(103415)

auc <- vector()

M <- 100

for(i in 1:M)
{
  for(snr in seq(0.1,1,0.01))
  {
    data <- try(with(parameters, simData(n    = n,
                                         coefficients = c[[1]],
                                         covariance   = 3.33,
                                         r.sqrd       = .1)))
    MAR_out <- try(with(parameters, makeMissing(data      = data,
                                                mechanism = "MAR",
                                                pm        = miss,
                                                preds     = colnames(data),
                                                snr       = snr
    )))
    auc <- c(auc, MAR_out$auc)
  }
}
df <- data.frame(cbind(auc, seq(0.1, 1, .01)))
colnames(df) <- c("auc", "snr")
df <- df[order(df$snr),]
plot(df)

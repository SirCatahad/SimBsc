## Clean up the workspace
rm(list=ls(all=TRUE))


## Libraries
library("mice")

source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(103415)

M <- 100
snr_seq <- seq(0.1, 1,0.01)

auc_snr_df <- data.frame(matrix(0,nrow=M*length(snr_seq), ncol=2))
colnames(auc_snr_df) <- c("auc", "snr")
auc_snr_df["snr"] <- snr_seq



cntr <- 1
for(i in 1:M)
{
  for(snr in snr_seq)
  {
    data <- try(with(parameters, simData(n    = n,
                                         coefficients = c[[1]],
                                         covariance   = 3.33,
                                         r.sqrd       = .1)))
    MAR_out_test <- try(with(parameters, makeMissing(data      = data,
                                                     mechanism = "MAR",
                                                     pm        = miss,
                                                     preds     = colnames(data),
                                                     snr       = snr
    )))
    auc_snr_df[cntr,"auc"] <- MAR_out_test$auc
    cntr <- cntr+1
  }
}
plot(auc_snr_df)
abline(v=.65,  col='red')
abline(h=.325, col='red')
abline(v=.75,  col='blue')
abline(h=.6,   col='blue')

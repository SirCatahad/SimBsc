## Clean up the workspace
rm(list=ls(all=TRUE))


## Libraries
library("mice")


source("init.R")
source("simMissingness.R")
source("sim_functions.R")


set.seed(103415)

test_iter <- 100
snr_seq <- seq(0.1, 1,0.01)

auc_snr_df <- data.frame(matrix(0,nrow=test_iter*length(snr_seq), ncol=2))
colnames(auc_snr_df) <- c("auc", "snr")
auc_snr_df["snr"] <- snr_seq



cntr <- 1
for(i in 1:test_iter)
{
  data <- try(with(parameters, simData(n    = n,
                                       beta   = 0)))
  for(snr in snr_seq)
  {
    
    MAR_out_test <- try(with(parameters, simLinearMissingness(pm       = miss,
                                                              data     = data,
                                                              snr      = snr,
                                                              type     = "high",
                                                              optimize = FALSE)
    ))
    auc_snr_df[cntr,"auc"] <- MAR_out_test$auc
    cntr <- cntr+1
  }
}
windowsFonts(A = windowsFont("Times New Roman"))
plot(auc_snr_df,
     family="A",
     xlab="AUC",
     ylab="SNR",cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
     yaxt='n',
     bty = "n",
     axes=FALSE)
axis(1, col="white", tcl=0)
axis(2, col="white", tcl=0, las=1)
abline(v=.65,  col='red')
abline(h=.325, col='red')
abline(v=.75,  col='blue')
abline(h=.6,   col='blue')

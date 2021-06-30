library("ggplot2")
library("gridExtra")

source("init.R")
source("sim_functions.R")


outputDir <- paste0("../data/study2/")

## Read in all the files
files <- list.files(path = outputDir, pattern = "\\.rds$", full.names = TRUE)
df <- lapply(files, readRDS)

## We now have a list of lists of data frames
df <- do.call("rbind", df)
##After this we have a list of data frames
df <- do.call("rbind", df)
##Now we have the combined data frame

## For some reason all saved data were characters, so we need to transform them
df[c("Slope X2", "conf.slope.X2.low", "conf.slope.X2.high", "Iteration")] <- sapply(df[c("Slope X2", "conf.slope.X2.low", "conf.slope.X2.high", "Iteration")], as.numeric)


## Check whether all iterations converge
all(is.na(df$Error))

##No errors found as the line above returned true, moving on

### Analysis and Visualization
## Bias estimations -----------------------------------------------------------------------------------------------------


bias_df <- df
## First, zero all conditions
complete_df <- data.frame(matrix(0,nrow=93,0))

## Bias is calculated as: 1/iter * sum_i=1_nsim ( theta_hat_i - theta ).
## First, we get each deviation
bias_df[grepl('cov3.33', bias_df$Condition),"Slope X2"] <- bias_df[grepl('cov3.33', bias_df$Condition),"Slope X2"] - betav[2]
bias_df[grepl('cov10', bias_df$Condition),"Slope X2"] <- bias_df[grepl('cov10', bias_df$Condition),"Slope X2"] - betav[3]

## Then we sum (per condition) and divide by nsum(iter)
complete_df[c("Condition", "Slope X2")] <- aggregate(bias_df["Slope X2"], by=bias_df["Condition"], FUN=(function(x) sum(x)/iter))

## MCSE(Bias) is calculated by: sqrt(1/nsim-1 * sum_i=1_nsim((theta_hat_i - theta_mean)^2)/nsim). This means we need a
## Theta_mean by condition
bias_df["Slope Mean Group"] <- ave(x=unlist(df["Slope X2"]), df["Condition"])


## Now calculate Monte carlo standard error of bias (MCSE(Bias))
complete_df[c("Condition", "Slope.SE")] <- aggregate(bias_df["Slope X2"], by=bias_df["Condition"],
                                                     FUN=function(x) sqrt((1/(iter-1)*sum((x^2)))/iter ))

## And finally Monte Carlo CIs
complete_df["slp.low"] <- sapply(complete_df["Slope X2"], FUN= (function(x,y) x - 1.96*y), complete_df["Slope.SE"])
complete_df["slp.high"] <- sapply(complete_df["Slope X2"], FUN= (function(x,y) x + 1.96*y), complete_df["Slope.SE"])

rm(bias_df)


## Coverage of confidence Intervals-----------------------------------------------------------------------------------------------------
## Coverage of confidence intervals: How often does the true population parameter lie in the
## estimated 95% confidence interval?

##New column that tells us whether the true population parameter lies in the confidence interval for each
##Iteration in each condition
df["CI_slp_cov"] <- NA

##Logical column whether the true parameter lies in the 95% confidence interval
df[grepl('cov0', df$Condition), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((betav[1] >= x && betav[1] <= y), 1, 0)),
  y=df[grepl('cov0', df$Condition), "conf.slope.X2.high"],
  x=df[grepl('cov0', df$Condition), "conf.slope.X2.low"])

df[grepl('cov3.33', df$Condition), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((betav[2] >= x && betav[2] <= y), 1, 0)),
  y=df[grepl('cov3.33', df$Condition), "conf.slope.X2.high"],
  x=df[grepl('cov3.33', df$Condition), "conf.slope.X2.low"])

df[grepl('cov10', df$Condition), "CI_slp_cov"] <- mapply(FUN=(function(x,y) 
  ifelse((betav[3] >= x && betav[3] <= y), 1, 0)),
  y=df[grepl('cov10', df$Condition), "conf.slope.X2.high"],
  x=df[grepl('cov10', df$Condition), "conf.slope.X2.low"])






##Percentage of CIs that contain the true population parameter
complete_df[c("Condition", "CI_slp_cov")] <- aggregate(df$CI_slp_cov, by=df["Condition"], FUN= function(x) (sum(x)))
complete_df["CI_slp_cov"] <- complete_df["CI_slp_cov"]/iter


complete_df["CI_slp_cov"] <- complete_df["CI_slp_cov"]*100


## SE of the slope CI coverage
complete_df["CI.slp.se"] <- as.vector(sapply(complete_df$CI_slp_cov, FUN=function(x) sqrt((x*(100-x))/iter)))

## Compute Monte Carlo CIs
complete_df["CIcov.slp.low"] <- as.vector(sapply(complete_df["CI_slp_cov"], FUN= (function(x,y) x - 1.96*y), complete_df$CI.slp.se ))
complete_df["CIcov.slp.high"] <- as.vector(sapply(complete_df["CI_slp_cov"], FUN= (function(x,y) x + 1.96*y), complete_df$CI.slp.se ))


## Empirical standard Error-----------------------------------------------------------------------------------------------------
## The standard deviation of β over 1,000 replications (henceforth the ‘empirical standard error’)

## Empirical SE is calculated as: sqrt(1/(nsim-1) * sum_i=1_nsim ( Theta_hat_i - Theta_mean)^2)
df["Slope_Mean"] <- ave(x=unlist(df["Slope X2"]), df["Condition"])

complete_df[c("Condition", "ESE")] <- aggregate(df["Slope X2"],
                                                by=df["Condition"],
                                                FUN=sd)



#Monte Carlo standard error of the empirical standard error
complete_df["MCSE_ESE"] <- as.vector(sapply(complete_df$ESE, FUN= function(x) x/(sqrt(2*(iter-1)))))

#complete_df[c("Condition", "MCSE_ESE")] <- aggregate(df$MCSE_ESE, by=df["Condition"], FUN= function(x) sum(x))

complete_df["SE.slp.low"] <- as.vector(sapply(complete_df["ESE"], FUN= (function(x,y) x - 1.96*y), complete_df$MCSE_ESE ))
complete_df["SE.slp.high"] <- as.vector(sapply(complete_df["ESE"], FUN= (function(x,y) x + 1.96*y), complete_df$MCSE_ESE ))

## Data prep-----------------------------------------------------------------------------------------------------

## Assign Colors to missing data conditions
complete_df["color"] <- NA

complete_df[grepl('complete_data', complete_df$Condition), "color"] <- 'Complete'
complete_df[grepl('MCAR', complete_df$Condition), "color"] <- 'MCAR'
complete_df[grepl('mcar', complete_df$Condition), "color"] <- 'MCAR'
complete_df[grepl('snr0.325', complete_df$Condition), "color"] <- 'Weak MAR'
complete_df[grepl('snr0.6', complete_df$Condition), "color"] <- 'Strong MAR'

complete_df["MatchType"] <- ""
complete_df[grepl('m1', complete_df$Condition), "MatchType"] <- "T1"
complete_df[grepl('m2', complete_df$Condition), "MatchType"] <- "T2"

complete_df["method"] <- NA
complete_df[grepl('complete_data', complete_df$Condition), "method"] <- "Complete"
complete_df[grepl('incomplete', complete_df$Condition), "method"] <- "Complete \nCases"
complete_df[grepl('pdmi', complete_df$Condition), "method"] <- "PDMI"
complete_df[grepl('pmm', complete_df$Condition), "method"] <- "PMM"

complete_df["k"] <- ""
complete_df[grepl('k1', complete_df$Condition), "k"] <- "k 1"
complete_df[grepl('k3', complete_df$Condition), "k"] <- "k 3"
complete_df[grepl('k5', complete_df$Condition), "k"] <- "k 5"
complete_df[grepl('k10', complete_df$Condition), "k"] <- "k10"

## Group properly
complete_df["Group"] <- complete_df["Condition"]
tmp <- read.csv("Grouping.csv")


complete_df$Group[match(tmp$Condition, complete_df$Condition)] <- tmp$X
complete_df$Group <- as.numeric(complete_df$Group)
complete_df <- complete_df[order(complete_df$Group),]

## Copy complete data cases twice for aesthetic purposes
complete_df <- rbind(complete_df, complete_df[grepl('complete_data', complete_df$Condition),], complete_df[grepl('complete_data', complete_df$Condition),])
complete_df <- complete_df[order(complete_df$Group),]

complete_df["Group"] <- 1:99



plot_df_a <- complete_df[grepl('cov0', complete_df$Condition),]
plot_df_b <- complete_df[grepl('cov3.33', complete_df$Condition),]
plot_df_c <- complete_df[grepl('cov10', complete_df$Condition),]

## Plots, drawings and everything neat -----------------------------------------------------------------------------------------------------


p11 <- makePlot(data    = plot_df_a, 
                xint    = 0, 
                title   = "ß = 0", 
                x       = plot_df_a$`Slope X2`, 
                xlow    = plot_df_a$slp.low, 
                xhigh   = plot_df_a$slp.high, 
                y       = plot_df_a$Group, 
                xlimits = c(-.2,.1))

p12 <- makePlot(data    = plot_df_b, 
                xint    = 0, 
                title   = "ß = 3.33", 
                x       = plot_df_b$`Slope X2`, 
                xlow    = plot_df_b$slp.low, 
                xhigh   = plot_df_b$slp.high, 
                y       = plot_df_b$Group, 
                xlimits = c(-1.5,1))

p13 <- makePlot(data    = plot_df_c, 
                xint    = 0, 
                title   = "ß = 10", 
                x       = plot_df_c$`Slope X2`, 
                xlow    = plot_df_c$slp.low, 
                xhigh   = plot_df_c$slp.high, 
                y       = plot_df_c$Group, 
                xlimits = c(-6,1.5))

grid_arrange_shared_legend(p11, p12, p13, nrow=1, ncol=3)

p21 <- makePlot(data    = plot_df_a, 
                xint    = 95, 
                title   = "ß = 0", 
                x       = plot_df_a$CI_slp_cov, 
                xlow    = plot_df_a$CIcov.slp.low, 
                xhigh   = plot_df_a$CIcov.slp.high, 
                y       = plot_df_a$Group, 
                xlimits = c(75,105))

p22 <- makePlot(data    = plot_df_b, 
                xint    = 95, 
                title   = "ß = 3.33", 
                x       = plot_df_b$CI_slp_cov, 
                xlow    = plot_df_b$CIcov.slp.low, 
                xhigh   = plot_df_b$CIcov.slp.high, 
                y       = plot_df_b$Group, 
                xlimits = c(40,100))

p23 <- makePlot(data    = plot_df_c, 
                xint    = 95, 
                title   = "ß = 10", 
                x       = plot_df_c$CI_slp_cov, 
                xlow    = plot_df_c$CIcov.slp.low, 
                xhigh   = plot_df_c$CIcov.slp.high, 
                y       = plot_df_c$Group, 
                xlimits = c(-5,100))


grid_arrange_shared_legend(p21, p22, p23, nrow=1, ncol=3)

p31 <- makePlot(data    = plot_df_a, 
                xint    = NULL, 
                title   = "ß = 0", 
                x       = plot_df_a$ESE, 
                xlow    = plot_df_a$SE.slp.low, 
                xhigh   = plot_df_a$SE.slp.high,  
                y       = plot_df_a$Group, 
                xlimits = c(.2,.55))

p32 <- makePlot(data    = plot_df_b, 
                xint    = NULL, 
                title   = "ß = 3.33", 
                x       = plot_df_b$ESE, 
                xlow    = plot_df_b$SE.slp.low, 
                xhigh   = plot_df_b$SE.slp.high,  
                y       = plot_df_b$Group, 
                xlimits = c(.25,.6))

p33 <- makePlot(data    = plot_df_c, 
                xint    = NULL, 
                title   = "ß = 10", 
                x       = plot_df_c$ESE, 
                xlow    = plot_df_c$SE.slp.low, 
                xhigh   = plot_df_c$SE.slp.high, 
                y       = plot_df_c$Group, 
                xlimits = c(.2,1.2))

grid_arrange_shared_legend(p31, p32, p33, nrow=1, ncol=3)

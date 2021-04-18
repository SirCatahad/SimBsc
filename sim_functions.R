##Libraries
library("SimDesign") ### KML: Are you actually using this library?

# Functions
## simData ------------------------------------------------------------------------------------------------------
simData <- function(n, coefficients, covariance, r.sqrd)
{

### KML: You only have 1 X variable. You can generate X as standard normal and
### get rid of all the marix operations. You don't need to worry about defining
### 'nrpred' or 'covariance'.
    
  #Construct covariance matrix
  nrpred <- length(coefficients) - 1
  sigma <- matrix(covariance, nrpred, nrpred)
  diag(sigma) <- 1.0
  
  #Generate data
  X <- rmvnorm(n = n, mean = rep(0, nrpred), sigma = sigma)

  #Generate error termn
  beta <- covariance ### KML: Why is beta = covariance? Where are you using
                     ### 'coefficients[2]'?
  if(r.sqrd > 0)
  {
    var.model <- t(beta) %*% cov(X) %*% beta
    var.residual <- (var.model/r.sqrd) - var.model
    U = rnorm(n, mean = 0, sd = sqrt(var.residual))
    #compute Y
    Y <- coefficients[1] + X  %*% beta + U
    data <- data.frame(X, Y)
  } else {
### KML: You can just simulate two independent normal variates here. The mean of
### Y is just the intercept and the variance is the residual variance.
      
    nrpred <- length(coefficients)
    sigma <- matrix(covariance, nrpred, nrpred)
    diag(sigma) <- 1.0
    data <- data.frame(rmvnorm(n = n, mean = rep(0, nrpred), sigma = sigma))
    colnames(data) <- c("X", "Y")
  }
  cov(data) ### KML: Don't put useless calculations in your functions.
  cov(X)
  #Return data
  data
}

## sim Data for quadratic function-------------------------------------------------------------------------------
simData2 <- function(n, coefficients, covariance, r.sqrd)
{
    ### KML: Same issue as above w.r.t. the single X variable. 
    
  #Construct covariance matrix
  nrpred <- length(coefficients) - 1
  sigma <- matrix(covariance, nrpred, nrpred)
  diag(sigma) <- 1.0
  
  #Generate data
  X <- rmvnorm(n = n, mean = rep(1, nrpred), sigma = sigma)
  
  #Generate error termn
  beta <- covariance
  if(r.sqrd > 0)
  {
### KML: You need to combine X and X^2 into a single predictor matrix before
### computing the residual variance
    var.model <- t(beta) %*% cov(X) %*% beta
    var.residual <- (var.model/r.sqrd) - var.model
    U = rnorm(n, mean = 0, sd = sqrt(var.residual))
                                        #compute Y

### KML: Where's your linear term?
    Y <- (X * X)  %*% beta + U
    data <- data.frame(X, Y)
  }
  else
  {
    nrpred <- length(coefficients)
    sigma <- matrix(covariance, nrpred, nrpred)
    diag(sigma) <- 1.0
    data <- data.frame(rmvnorm(n = n, mean = rep(1, nrpred), sigma = sigma))
    colnames(data) <- c("X", "Y")
    #data["Y"] <- 2*(data["Y"]*data["Y"])
  }
  
  
  #Return data
  data
}

## make Missing -----------------------------------------------------------------------------------------------------
# data       - the data frame which should get missing observations
# mechanism  - the mechanism of missing data, by default MCAR
# percent    - the proportion of observations that should be set to missing (NA)
# indices    - A vector of indices indicating which columns should contain missing values

### KML: This function doesn't make much sense. You're returning the response
### vector, which is the same thing returned by the simLinearMissingness()
### function. So, why what simLinearMissingness() in another function with the
### same return value? You need a function to generate MCAR missingness, but it
### doesn't have to also generate MAR missingness. If you create a combined
### function, it should take some complete data and return incomplete data.


makeMissing <- function(data, 
                        mechanism="MCAR", 
                        pm, 
                        preds, 
                        snr=NULL)
{
  #MAR missing data mechanism
  if(mechanism=="MAR")
  {
    out <- simLinearMissingness(pm       = pm,
                                data     = data,
                                snr      = snr,
                                preds    = preds,
                                type     = "high",
                                optimize = FALSE)
    
    out
  }
  
  #MCAR missing data mechanism
  else if(mechanism=="MCAR")
  {
    r <- sample(1:nrow(data), nrow(data)*pm)
    tmp <- rep(FALSE, nrow(data))
    tmp[r] <- TRUE
    r <- tmp
    out <- list(r   = r)#,
                #eta = eta2,
                #auc = auc,
                #snr = sd(eta) / sqrt(var(eta2) - var(eta)))
    #return
    out
  }
  else
  {
    stop("Undefined or unsupported missing data mechanism.")
  }
}
## zip function -------------------------------------------------------------------------------------------------------------------
## R implementation of Pythons zip() function
zip <- function(...) 
{
  mapply(list, ..., SIMPLIFY = FALSE)
}

## Analyze (probably unnecessary)------------------------------------------------------------------------------------------------
## Analyze
## Extracts relevant information from datasets, namely:
## - Intercept and slope coefficients
## - confidence intervals
## - Standard errors

### KML: This function isn't really fit for purpose, for two reasons:
###      1. It won't work for MI data. You can use the with.mids() function to
###         fit the models to the imputed data.
###      2. You're saving a bunch of stuff that you don't need. You only need
###         the coefficients and CIs.

analyze <- function(data)
{
  #Remove NAs
  data <- data[complete.cases(data),]
  #run regression model
  reg <- lm(Y ~ ., data = data)
  
  #build output
  out <- list(
    cov           = cov(data),
    coefficients  = reg$coefficients,
    residuals     = reg$residuals,
    regsum        = summary(reg),
    c_int         = confint(reg)
  )
  
  out
}

### KML: I haven't checked any of the analysis/visualization code yet. You need
### to get the simulation running correctly before worrying about visualizing
### the results.

makePlot <- function(data, xint, title, x, xlow, xhigh, y, xlimits)
{
  plt <- ggplot(data, aes(col=color,group=method, xmin=xlow, x=x, xmax=xhigh))+
    geom_pointrange(aes(y = y))+
    #geom_linerange(aes(x=Intercept, ymin=Condition, ymax=Condition, xmin=conf.intercept.low, xmax=conf.intercept.high)) +
    labs(title=title)+
    #  coord_flip()+
    theme(plot.title = element_text(hjust = 0.5),
          strip.text.y.left = element_text(angle = 0),
          panel.spacing.y=unit(.5,"lines"),
          # Remove background
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #panel.background = element_blank(), axis.line = element_line(colour = "black"),
          #Remove y ticks and labels
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          #Place labels outside
          strip.placement = "outside",
          strip.background =element_rect(fill="white"),
          #Place legend at the top
          legend.position="top"
    )+
    scale_color_manual(values=c("black","red","green","blue")) +
    scale_x_continuous(limits=xlimits)+
    facet_grid(rows=vars(method,k),
               scales="free",
               space="free_y",
               drop=TRUE,
               switch="y",
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 list(do.call(paste, c(labels, list(sep = "\n"))))
               })+ 
    geom_vline(xintercept = xint)
  plt
}

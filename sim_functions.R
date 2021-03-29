##Libraries
library("SimDesign")

# Functions
##-------------------------------------------------------------------------------------------------------------------##
simData <- function(n, coefficients, covariance, r.sqrd)
{
  #Construct covariance matrix
  nrpred <- length(coefficients) - 1
  sigma <- matrix(covariance, nrpred, nrpred)
  diag(sigma) <- 1.0
  
  #Generate data
  X <- rmvnorm(n = n, mean = rep(0, nrpred), sigma = sigma)

  #Generate error termn
  beta <- coefficients[-1]
  if(r.sqrd > 0)
  {
    var.model <- t(beta) %*% cov(X) %*% beta
    var.residual <- (var.model/r.sqrd) - var.model
    U = rnorm(n, mean = 0, sd = sqrt(var.residual))
    #compute Y
    Y <- coefficients[1] + X  %*% beta + U
    data <- data.frame(X, Y)
  }
  else
  {
    nrpred <- length(coefficients)
    sigma <- matrix(covariance, nrpred, nrpred)
    diag(sigma) <- 1.0
    data <- data.frame(rmvnorm(n = n, mean = rep(0, nrpred), sigma = sigma))
    colnames(data) <- c("X", "Y")
  }
  
  
  #Return data
  data
}

##-------------------------------------------------------------------------------------------------------------------##
simData2 <- function(n, coefficients, covariance, r.sqrd)
{
  #Construct covariance matrix
  nrpred <- length(coefficients) - 1
  sigma <- matrix(covariance, nrpred, nrpred)
  diag(sigma) <- 1.0
  
  #Generate data
  X <- rmvnorm(n = n, mean = rep(1, nrpred), sigma = sigma)
  
  #Generate error termn
  beta <- coefficients[-1]
  if(r.sqrd > 0)
  {
    var.model <- t(beta) %*% cov(X) %*% beta
    var.residual <- (var.model/r.sqrd) - var.model
    U = rnorm(n, mean = 0, sd = sqrt(var.residual))
    #compute Y
    Y <- coefficients[1] + X*X  %*% beta + U
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

##-------------------------------------------------------------------------------------------------------------------##
# data       - the data frame which should get missing observations
# mechanism  - the mechanism of missing data, by default MCAR
# percent    - the proportion of observations that should be set to missing (NA)
# indices    - A vector of indices indicating which columns should contain missing values
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
##-------------------------------------------------------------------------------------------------------------------##
## R implementation of Pythons zip() function
zip <- function(...) 
{
  mapply(list, ..., SIMPLIFY = FALSE)
}

##-------------------------------------------------------------------------------------------------------------------##
## Analyze
## Extracts relevant information from datasets, namely:
## - Intercept and slope coefficients
## - confidence intervals
## - Standard errors
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

makePlot <- function(data, xint, title, x, xlow, xhigh, y, xlimits)
{
  plt <- ggplot(data, aes(col=color,group=method, xmin=xlow, x=x, xmax=xhigh))+
    geom_pointrange(aes(y = y))+
    #geom_linerange(aes(x=Intercept, ymin=Condition, ymax=Condition, xmin=conf.intercept.low, xmax=conf.intercept.high)) +
    labs(title=title)+
    #  coord_flip()+
    theme(plot.title = element_text(hjust = 0.5),
          strip.text.y.left = element_text(angle = 0),
          panel.spacing.y=unit(0,"lines"),
          # Remove background
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
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
               scales="fixed",
               drop=TRUE,
               switch="y",
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 list(do.call(paste, c(labels, list(sep = "\n"))))
               })+ 
    geom_vline(xintercept = xint)
  plt
}

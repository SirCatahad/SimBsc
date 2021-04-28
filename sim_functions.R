# Functions
## simData ------------------------------------------------------------------------------------------------------
simData <- function(n, covariance, mod=NULL)
{
  
  #In the second study, the data generation model equals beta*x^2
  #In the first it is beta*x.
  if(is.null(mod))
  {
    X <- data.frame(rnorm(n = n, mean = 0, sd = 1))
    
  } else {
    X <- data.frame(rnorm(n = n, mean = 1, sd = 1))
  }
  
  #Beta is equal to the covariance
  beta <- covariance 
  
  #Generate error termn
  U = rnorm(n, mean=0, sd=sqrt(100))
  
  #compute Y, Intercept will be zero.
  if(is.null(mod))
  {
    Y <- beta * X + U
  } else {
    
    Y <-   beta * X^2 + U
  }
  
  data <- data.frame(X, Y)
  colnames(data) <- c("X", "Y")

  ## Return
  data
}

## make Missing -----------------------------------------------------------------------------------------------------
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
    tmp <- simLinearMissingness(pm       = pm,
                                data     = data,
                                snr      = snr,
                                preds    = preds,
                                type     = "high",
                                optimize = FALSE)
    
    data[tmp$r,1] <- NA
    data
    
  }
  
  #MCAR missing data mechanism
  else if(mechanism=="MCAR")
  {
    r <- sample(1:nrow(data), nrow(data)*pm)
    tmp <- rep(FALSE, nrow(data))
    tmp[r] <- TRUE
    data[tmp,1] <- NA
    data
  }
  else
  {
    stop("Undefined or unsupported missing data mechanism.")
  }
}

## Analyze------------------------------------------------------------------------------------------------

analyze <- function(data, study)
{
  #Remove NAs
  data <- data[complete.cases(data),]
  #run regression model
  if(study == "study1")
  {
    reg <- lm(Y ~ X, data = data)
  } else {
    reg <- lm(Y ~ X + I(X^2), data = data)
  }
  
  
  #build output
  out <- c(reg$coefficients[-1],
           confint(reg)[-1,c(1,2)])
  out
}

analyze_MI <- function(MI_object, study)
{
  if(study=="study1")
  {
    fit <- with(MI_object, lm(Y ~ X))
  } else {
    fit <- with(MI_object, lm(Y ~ X + I(X^2)))
  }
  
  
  result <- summary(pool(fit), conf.int = TRUE, conf.level = .95)
  
  #build output
  out <- c(result$estimate[-1],
           result$`2.5 %`[-1],
           result$`97.5 %`[-1])
  out
}


##
## Visualization stuff ---------------------------

### KML: I haven't checked any of the analysis/visualization code yet. You need
### to get the simulation running correctly before worrying about visualizing
### the results.

### BP : True, on it

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

## Run one iteration --------------------------------

runIteration <- function(covariances, parameters, km, mtype, snr, study, iter)
{
  
  #Data generation and storage
  if(study=="study1")
  {
    storage <- data.frame(matrix(0,
                                 nrow= 93, #nr of conditions.
                                 ncol=  5  #saved variables
    ))
    colnames(storage) <- c("Slope X",
                           "conf.slope.X.low",
                           "conf.slope.X.high",
                           "Condition",
                           "Iteration")
  } else {
    storage <- data.frame(matrix(0,
                                 nrow= 93, #nr of conditions.
                                 ncol=  8  #saved variables
    ))
    colnames(storage) <- c("Slope X",
                           "Slope X2",
                           "conf.slope.X.low",
                           "conf.slope.X2.low",
                           "conf.slope.X.high",
                           "conf.slope.X2.high",
                           "Condition",
                           "Iteration")
  }
  
  
  ## Storing all important variables plus an identifier so we can easily order by conditions later on
  
  con <- 1
  
  
  
  for(cv in covariances) 
  {
    if(study=="study1")
    {
      #Generate data with given parameters
      data <- try(with(parameters, simData(n    = n,
                                           covariance   = cv)))
    } else {
      #Generate data with given parameters
      data <- try(with(parameters, simData(n    = n,
                                           covariance   = cv,
                                           "Quadratic")))
    }
    
    ## Get relevant information
    storage[con, ] <- c(analyze(data, study), paste0("complete_data_cov",cv), iter)
    ## Update counter
    con <- con+1 

    #Poke holes into the dataset
    data_MCAR <- try(with(parameters, makeMissing(data      = data,
                                                 mechanism = "MCAR",
                                                 pm        = miss
    )))
    
    #Store
    storage[con, ] <- c(analyze(data_MCAR, study), paste0("incomplete_MCAR_cov",cv), iter)
    con <- con+1
    
    PDMI <- try(mice(data   = data_MCAR,
                     m      = 10,
                     method = "norm",
                     printFlag = FALSE
    ))
    #Storage
    storage[con, ] <- c(analyze_MI(PDMI, study), paste0("pdmi_mcar_cov",cv), iter)
    
    con <- con+1
    
    conds <- expand.grid(mtype, km)
    
    ## Impute with PMM
    for(cs in 1:nrow(conds))
    {
      PMM <- try(mice(data   = data_MCAR,
                      m      = 10,
                      method = "pmm",
                      matchtype = conds[cs,1],
                      donors = conds[cs,2],
                      printFlag = FALSE
      ))

      
      
      #Store
      storage[con, ] <- c(analyze_MI(PMM, study), paste0("pmm_mcar_cov",cv, "_k", conds[cs,2],"_m", conds[cs,1]), iter)
      con <- con+1
    }
    
    for(s in snr)
    {
      data_MAR <- try(with(parameters, makeMissing(data      = data,
                                                  mechanism = "MAR",
                                                  pm        = miss,
                                                  preds     = colnames(data),
                                                  snr       = s
      )))
      
      #Store
      storage[con, ] <- c(analyze(data_MAR, study), paste0("incomplete_mar_cov",cv,"_snr",s), iter)
      con <- con+1
      
      PDMI <- try(mice(data   = data_MAR,
                       m      = 10,
                       method = "norm",
                       printFlag = FALSE
      ))
      #Storage
      storage[con, ] <- c(analyze_MI(PDMI, study), paste0("pdmi_mar_cov",cv,"_snr",s), iter)
      
      con <- con+1
      
      ## Impute with PMM
      for(cs in 1:nrow(conds))
      {
        PMM <- try(mice(data   = data_MAR,
                        m      = 10,
                        method = "pmm",
                        matchtype = conds[cs,1],
                        donors = conds[cs,2],
                        printFlag = FALSE
        ))

        #Store
        storage[con, ] <- c(analyze_MI(PMM, study), paste0("pmm_mar_cov",cv,"_snr",s,"_k", conds[cs,2],"_m", conds[cs,1]), iter)
        con <- con+1
      }
    }
  }
  storage
}


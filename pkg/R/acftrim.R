#####################
# acf using x trimming
# input
# x: time series without NA as vector
# lag.max: the maximal lag of interest
# trim: the part of largest and smallest observations, which should be trimmed
# output: autocorrelation function
#####################

acftrim <- function(x, lag.max, trim=0.1) {
  n <- length(x)
  lags <- 1:lag.max
  
  # protective checks:
  if (trim < 0) {
  	warning("The trimming proportion has to be positive and is therefore set to 0.1 instead")
  	trim <- 0.1
  	}
  if (trim > 0.5) {
  	warning("You can not trim more then the half smalles and largest values. The trimming constant is therefore set to 0.1.")
  	trim <- 0.1
  	}
  
  # trimming the time series:  
  qu <- quantile(x, trim)
  qo <- quantile(x, 1-trim)  
  Lu <- x>= qu
  Lo <- x<= qo
  L <- Lo*Lu			# untrimmed observations 
  if (sum(L) < 1) {
  	stop("There are no untrimmed observations left to calculate the acf")
 	}
    
  # calculation acf:
  Xquer <- sum(x*L)/sum(L)	# time series trimmed mean
  gamma <- numeric(length(lags))
  for (i in lags) {
    X1 <- x[1:(n-i)]
    L1 <- L[1:(n-i)]
    X2 <- x[(i+1):n]
    L2 <- L[(i+1):n]
    gamma[i] <- sum((X1-Xquer)*(X2-Xquer)*L1*L2)/(sum(L1*L2))
  }  
  if (sum(L1*L2)==0) {
  	warning("It have been trimmed too many observations, such that it was not possible to calculate the acf for all lags")
 	}
  gamma0 <- sum((x-Xquer)^2*L)/sum(L)	# time series trimmed variance
  acfvalues_biased <- gamma/gamma0
  
  # transformation for unbiasedness:
  load(system.file("extdata", "trimsimv", package = "robts"))	# loading the simulated expected values
  acfvalues <- sapply(acfvalues_biased, linearinterpol, a=expectations, b=values)
  
  return(acfvalues)
}

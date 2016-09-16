#####################
# calculates the acf using a RA-approach
# input
# x: time series without NA as vector
# lag.max: maximal lag of interest
# Psi: which Psi-Function should be used? (Huber and Tukey are available)
# locfn: a mean estimator
# scalefn: an estimator of scatter 
# ...: Robustness-Parameter for Psi-Functions
# output: acf
#####################

acfRA <- function(x, lag.max, Psi="Huber", locfn=median, scalefn=mad, ...) {
  n <- length(x)
  lags <- 1:lag.max
  
  # protective measures:
  if (!is.function(locfn)) {
  	warning("This is no suitable location estimator. Median is used instead")
  	locfn <- median
 	}
  if (!is.function(scalefn)) {
  	warning("This is no suitable scatter estimator. MAD is used instead")
  	scalefn <- mad
 	}
  
  # calculation of mean and variance:
  locval <- try(locfn(x), silent=TRUE)
  if (inherits(locval, "try-error")) {
  	warning("Something went wrong with the mean estimation")
  	return(NA)
 	}
  scaleval <- scalefn(x)
  if (inherits(scaleval, "try-error")) {
  	warning("Something went wrong with the scale estimation")
  	return(NA)
 	}
  
  # Transformation of values:
  x <- (x-locval)/scaleval
  if (Psi=="Huber") {
  	x <- Huber(x,...)
 	}
  if (Psi=="Tukey") {
  	x <- tukeypsi(x,...)
 	}
  if (!any(Psi==c("Huber","Tukey")))
  	{warning("This is no suitable Psi-function. Huber is used instead.")
  	x <- Huber(x)
 	}
  
  # calculation of acf (biased!): 
  acfvalues_biased <- acf(x, demean=FALSE, plot=FALSE, lag.max=lag.max)$acf[-1]

  # transformation for unbiasedness:
  if (Psi=="Huber") load(system.file("extdata", "rahusimv", package = "robts")) # loading the simulated expected values
  if (Psi=="Tukey")	load(system.file("extdata", "ratusimv", package = "robts")) # loading the simulated expected values
  if (!any(Psi==c("Huber", "Tukey"))) {
    acfvalues <- sapply(acfvalues_biased, linearinterpol, a=expectations, b=values)
  	warning("The estimation could be slightly biased")
 	}
 	
  return(acfvalues)
}
#####################
# calculates the acf using a RA-approach
# input
# x: time series without NA as vector
# lag.max: maximal lag of interest
# psi: which Psi-Function should be used? ("huber" and "bisquare" are available)
# locfn: a mean estimator
# scalefn: an estimator of scatter 
# ...: Robustness-Parameter for Psi-Functions
# output: acf
#####################

acfRA <- function(x, lag.max, psi="huber", locfn=median, scalefn=mad, ...) {
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
  x_transformed <- M_psi((x-locval)/scaleval, psi=psi, ...)
  # calculation of acf (biased!): 
  acfvalues_biased <- acf(x_transformed, demean=FALSE, plot=FALSE, lag.max=lag.max)$acf[-1]

  # transformation for unbiasedness:
  if (psi=="huber") load(system.file("extdata", "rahusimv", package = "robts")) # loading the simulated expected values
  if (psi=="bisquare") load(system.file("extdata", "ratusimv", package = "robts")) # loading the simulated expected values
  if (psi %in% c("huber", "bisquare")) {
    acfvalues <- sapply(acfvalues_biased, linearinterpol, a=expectations, b=values)
  }else{
    warning("No bias correction has been made. The estimation could be slightly biased")
 	}
 	
  return(acfvalues)
}
###################
# calculates the acf based on median correlation
# input
# x: time series without NAs as vector
# lag.max: maximal lag of interest
# output: calculated acf
###################


acfmedian <- function(x, lag.max) {
  n <- length(x)
  lags <- 1:lag.max
  
  # calculating the acf (biased!):
  x_centered <- x - median(x)
  acfvalues_biased <- numeric(length(lags))
  for (i in lags) {
  	acfvalues_biased[i] <- mediancor(x_centered[1:(n-i)], x_centered[(i+1):n])
  }
  
  # transformation for unbiasedness:
  load(system.file("extdata", "chack2", package = "robts")) # loading the simulated expected values for the median correlation
 	acfvalues <- sapply(acfvalues_biased, linearinterpol, a=expectations, b=values)
 	
  return(acfvalues)
}
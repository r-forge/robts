###################
# calculates the acf using the GK approach
# input
# x: time series (without NA) as vector
# lag.max: the maximal lag of interest
# scalefn: which variance estimator function to use (Qn, mad, scaleTau2, reweightedQn, ...)
# ...: passed to function scalefn
# output: autocorrelation function
###################

acfrob.GK <- function(x, lag.max, scalefn = Qn, ...) {
  n <- length(x)
  lags <- 1:lag.max
  
  # calculating the acf:
  acfvalues <- numeric(length(lags))
  for (i in lags) {
    acfvalues[i] <- corGK(x[1:(n-i)], x[(i+1):n], scalefn=scalefn, ...)
  }
 
  return(acfvalues)
}

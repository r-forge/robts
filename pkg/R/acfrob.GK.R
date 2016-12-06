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
  
  are <- NA #factor is NA unless the scale function matches one of the following alternatives:
  if(identical(scalefn, Qn)) are <- sqrt(1/0.8227)
  if(identical(scalefn, scaleTau2)) are <- sqrt(1/0.8)
  if(identical(scalefn, mad)) are <- sqrt(1/0.3674)
  if(identical(scalefn, sd)) are <- 1

 	res <- list(
   acfvalues = acfvalues,
   are = are
  )
  	
  return(res)
}

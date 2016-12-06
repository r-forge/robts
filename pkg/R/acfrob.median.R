###################
# calculates the acf based on median correlation
# input
# x: time series without NAs as vector
# lag.max: maximal lag of interest
# output: calculated acf
###################


acfrob.median <- function(x, lag.max, biascorr = TRUE) {
  n <- length(x)
  lags <- 1:lag.max
  
  # calculating the acf (biased!):
  x_centered <- x - median(x)
  mediancor <- function(x, y) median(x*y)/median(x^2) #median correlation
  acfvalues_biased <- numeric(length(lags))
  for (i in lags) {
  	acfvalues_biased[i] <- mediancor(x_centered[1:(n-i)], x_centered[(i+1):n])
  }
  if(!biascorr) return(acfvalues_biased)
  
  # transformation for unbiasedness:
  load(system.file("extdata", "acfbiascorr_median", package = "robts")) # loading the simulated expected values for the median correlation
 	acfvalues <- sapply(acfvalues_biased, linearinterpol, a=get("expectations"), b=get("values"))
 	
 	res <- list(
   acfvalues = acfvalues,
   are = NA
  )
  	
  return(res)
}
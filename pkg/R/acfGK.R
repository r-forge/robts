###################
# calculates the acf using the GK approach
# input
# timeseries: timeseries (without NA) as vector
# maxlag: the maximal lag of interest
# method: which variance estimator to use (Qn, MAD and Tau available)
# output: autocorrelation function
###################

acfGK <- function(timeseries,maxlag,GK.method="Qn",...) {

n <- length(timeseries)

# calculating the acf

acfvalues <- numeric(maxlag)
for (i in 1:maxlag) {
acfvalues[i] <- GK(timeseries[1:(n-i)],timeseries[(i+1):n],method=GK.method,...)
}
return(acfvalues)
}

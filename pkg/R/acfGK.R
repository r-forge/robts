###################
# calculates the acf using the GK approach
# input
# timeseries: timeseries (without NA) as vector
# maxlag: the maximal lag of interest
# method: which variance estimator to use (Qn, MAD and Tau available)
# output: autocorrelation function
###################

acfGK <- function(timeseries,maxlag,method=c("Qn", "MAD", "Tau")) {

# protective measures

if(sum(is.na(timeseries))>0) {
	warning("There are NA in your timeseries you should use a procedure to replace this values first.")
	return(NA)
	}
n <- length(timeseries)
if (n< 2+maxlag) {
	warning("You can't calculate this amount of lags. We will compute the maximal possible lag of n-2")
	maxlag <- n-2
	}
if (n< 4*maxlag) {
	warning("Brockwell and Davis suggest to calculate only lags less n/4. Nevertheless we will calculate all lags you want, but you should be aware that the estimated acf for higher lags could be unreasonable.")
	}

# calculating the acf

acfvalues <- numeric(maxlag)
for (i in 1:maxlag) {
acfvalues[i] <- GK(timeseries[1:(n-i)],timeseries[(i+1):n],method=method)
}
return(acfvalues)
}
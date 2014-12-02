###################
# calculates the acf based on mediancorrelation
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximal lag of interest
# output: calculated acf
###################


acfmedian <- function(timeseries,maxlag) {

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

# calculating the acf (biased!)

timeseries <- timeseries - median(timeseries)
acfvalues <- numeric(maxlag)
for (i in 1:maxlag) {
	acfvalues[i] <- mediancor(timeseries[1:(n-i)],timeseries[(i+1):n])
}

# transformation for unbiasedness

AA <- get(load(system.file("extdata", "chack2", package = "robts")))				# loading the simulated expection values for the mediancorrelation
b <- seq(from=-1,to=1,by=0.01)		# lattice where the mediancorrelation was simulated
acfvalues2 <- numeric(maxlag)
for (i in 1:maxlag) {
	acfvalues2[i] <- linearinterpol(acfvalues[i],AA,b)
	}
return(acfvalues2)
}
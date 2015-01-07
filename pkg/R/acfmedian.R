###################
# calculates the acf based on mediancorrelation
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximal lag of interest
# output: calculated acf
###################


acfmedian <- function(timeseries,maxlag) {

n <- length(timeseries)


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
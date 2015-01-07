#####################
# acf using timeseries trimming
# input
# timeseries: timeseries without NA as vector
# maxlag: the maximal lag of interest
# alpha: the part of largest and smallest observations, which should be trimmed
# output: autocorrelation function
#####################

acftrim <- function(timeseries,maxlag,alpha=0.1) {

# protective measures

n <- length(timeseries)

if (alpha < 0) {
	warning("The trimming Quantilke has to be positive and is therefore set to 0.1 instead")
	alpha <- 0.1
	}
if (alpha > 0.5) {
	warning("You can not trim more then the half smalles and largest values. The trimming constant is therefore set to 0.1.")
	alpha <- 0.1
	}

# trimming the timeseries

qu <- quantile(timeseries,alpha)
qo <- quantile(timeseries,1-alpha)

Lu <- timeseries>= qu
Lo <- timeseries<= qo
L <- Lo*Lu			# untrimmed observations

if (sum(L) < 1) {
	warning("There are no untrimmed observations left to calculate the acf")
	return(NA)
	}

Xquer <- sum(timeseries*L)/sum(L)	# timeseries trimmed mean

# calculation acf

gamma <- numeric(maxlag)
for (i in 1:maxlag) {
X1 <- timeseries[1:(n-i)]
L1 <- L[1:(n-i)]
X2 <- timeseries[(i+1):n]
L2 <- L[(i+1):n]
gamma[i] <- sum((X1-Xquer)*(X2-Xquer)*L1*L2)/(sum(L1*L2))
}

if (sum(L1*L2)==0) {
	warning("You have trimmed so many observations, that it was not possible to calculate the acf for all lags")
	}

gamma0 <- sum((timeseries-Xquer)^2*L)/sum(L)	# time series trimmed variance

acfvalues <- gamma/gamma0

AA <- get(load(system.file("extdata", "trimsimv", package = "robts")))				# loading the simulated expection values for the mediancorrelation
b <- seq(from=-1,to=1,by=0.01)		# lattice where the mediancorrelation was simulated
acfvalues2 <- numeric(maxlag)
for (i in 1:maxlag) {
	acfvalues2[i] <- linearinterpol(acfvalues[i],AA,b)
	}
return(acfvalues2)
}

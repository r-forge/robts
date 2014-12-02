#####################
# calculates the acf using a RA-approach
# input
# timeseries: timeseries without NA as vector
# maxlag: maximal lag of interest
# Psi: which Psi-Function should be used? (Huber and Tukey are available)
# meanvalue: a mean estimator
# scattervalue: an estimator of scatter 
# ...: Robustness-Parameter for Psi-Functions
# output: acf
#####################

acfRA <- function(timeseries,maxlag,Psi="Huber",meanvalue=median,scattervalue=mad,...) {

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
if (!is.function(meanvalue)) {
	warning("This is no suitable location estimator. Median is used instead")
	meanvalue <- median
	}
if (!is.function(scattervalue)) {
	warning("This is no suitable scatter estimator. MAD is used instead")
	scattervalue <- mad
	}

# calculation of mean and variance
meanv <- try(meanvalue(timeseries),silent=TRUE)
if (inherits(meanv,"try-error")) {
	warning("Something went wrong with the mean estimation")
	return(NA)
	}
scatterv <- scattervalue(timeseries)
if (inherits(scatterv,"try-error")) {
	warning("Something went wrong with the mean estimation")
	return(NA)
	}

# Transformation of values

timeseries <- (timeseries-meanv)/scatterv
if (Psi=="Huber") {
	timeseries <- Huber(timeseries,...)
	}
if (Psi=="Tukey") {
	timeseries <- tukeypsi(timeseries,...)
	}
if (!any(Psi==c("Huber","Tukey")))
	{warning("This is no suitable Psi-function. Huber is used instead.")
	timeseries <- Huber(timeseries)
	}

# calculation of acf

acfvalues <- acf(timeseries,demean=FALSE,plot=FALSE,lag.max=maxlag)$acf[-1]
if (Psi=="Huber") {
	AA <- get(load(system.file("extdata", "rahusimv", package = "robts")))		# loading the simulated expection values for the mediancorrelation
	b <- seq(from=-1,to=1,by=0.01)		# lattice where the mediancorrelation was simulated
	acfvalues2 <- numeric(maxlag)
	for (i in 1:maxlag) {
	acfvalues2[i] <- linearinterpol(acfvalues[i],AA,b)
	}
	}
if (Psi=="Tukey") {
	AA <- get(load(system.file("extdata", "ratusimv", package = "robts")))				# loading the simulated expection values for the mediancorrelation
	b <- seq(from=-1,to=1,by=0.01)		# lattice where the mediancorrelation was simulated
	acfvalues2 <- numeric(maxlag)
	for (i in 1:maxlag) {
	acfvalues2[i] <- linearinterpol(acfvalues[i],AA,b)
	}
	}
if (!any(Psi==c("Huber","Tukey"))) {
	acfvalues2 <- acfvalues
	warning("The estimation could be slightly biased")
	}
return(acfvalues2)
}
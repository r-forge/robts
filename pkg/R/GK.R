###############
# auxiliary function: calculates (2-dim) correlation based on variances (GK-approach)
# input
# x: first variable to calculate correlation
# y: second variable to calculate correlation
# ...: parameters for Tau estimator
# method: what variance estimator to use (Qn, Tau and MAD available)
# output: correlation-coefficient
###############
varS <- function(x,boundq=0.999) {
n <- length(x)
s0 <- Qn(x)
m0 <- median(x)
robmalahanobis <- ((x-m0)/s0)^2
ordermalahanobis <- sort(robmalahanobis)

bound <- qchisq(boundq,1)			# bound from which distribution function of malahanobis distances is compared
dn <- max((pchisq(ordermalahanobis,df=1)-rank(ordermalahanobis)/n+1/n)*(ordermalahanobis>=bound)) # maximum positiv! difference between theoretical and empirical distribution function, compared from borderquantile (+1/n is result of jumpdiscontinuity of empirical distribution) 

alphan <- n-floor(dn*n)		# calculation of the cutoff value
alphan <- ordermalahanobis[alphan]
weight <- (robmalahanobis < alphan)	 # cut all observations which are greater  then the observation where the difference between the distribution functions was maximal
gooddata <- x[weight==1]
return(sd(gooddata))
}

GK <- function(x,y,method="Qn",...) {

### choose varaince estimator

if (!any(method==c("Qn","Tau","MAD","effi"))) 
	{warning("This is no suitable variance estimator. Qn is used instead.")
	varest <- Qn}
if (method=="Qn") {
	varest <- Qn
	}
if (method=="Tau") {
	varest <- function(x) scaleTau2(x,...)
	}
if (method=="MAD") {
	varest <- mad
	}
if (method=="effi"){
	varest <- varS
	}

### calculating necessary variances

varsum <- varest(x+y)^2
vardif <- varest(x-y)^2
if(varsum+vardif==0){warning("Something is wrong with your timeseries since variance estimation is 0.")
	return(NA)	
	}
return((varsum-vardif)/(varsum+vardif))
}

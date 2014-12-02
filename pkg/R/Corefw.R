#####################
# auxiliary function: calculates multivariate correlation (a fully efficient reweighting step)
# input
# data: multivariate data as matrix (observations in rows)
# boundq: Quantile-bound for comparism of theoretical and empirical distribution function
# startestimator: which startestimator should be used (raw MCD, reweighted MCD and S estimator available)
# output: correlation matrix
######################


Corefw <- function(data,boundq=0.975,startestimator="wMCD") {

# dimensions

n <- length(data[,1])
p <- length(data[1,])

# protective measureas

if (boundq<=0) {
warning("Arguments for the quantilefunction have to be positiv. A 0.975-quantile is used instead.")
boundq <- 0.975
}

if (boundq< 0.5) {
warning("Using a quantile smaller than the median is usually not reasonable. Nonetheless we will use it.")
}

if (boundq> 1) {
warning("Arguments for the quantilefunction have to be smaller then 1. A 0.975-quantile is used instead.")
boundq <- 0.975
}

# start estimation

if (startestimator=="wMCD") {
	estimate <- try(covMcd(data),silent=TRUE)
	if (inherits(estimate,"try-error")) {
		warning("Calculation of MCD failed")
		return(NA)
		}
	sigma <- estimate$cov
	center <- estimate$center
	}
if (startestimator=="rMCD") {
	estimate <- try(covMcd(data),silent=TRUE)
	if (inherits(estimate,"try-error")) {
		warning("Calculation of MCD failed")
		return(NA)
		}
	sigma <- estimate$raw.cov
	center <- estimate$raw.center 
	}
if (startestimator=="S") {
	estimate <- try(CovSest(data),silent=TRUE)
	if (inherits(estimate,"try-error")) {
		warning("Calculation of S estimator failed")
		return(NA)
		}
	sigma <- estimate@cov
	center <- estimate@center
	}
if (!any(startestimator==c("wMCD","rMCD","S"))) {
	warning("This is no suitable start estimator. We will use the weighted MCD instead.")
	estimate <- try(covMcd(data),silent=TRUE)
	if (inherits(estimate,"try-error")) {
		warning("Calculation of MCD failed")
		return(NA)
		}
	sigma <- estimate$cov
	center <- estimate$center
	}

try(inver <- solve(sigma),silent=TRUE)
if (inherits(inver,"try-error")) {
	warning("Correlationmatrix is not invertable. Therfore fully efficient reweighting is not possible, return start estimation instaed.")
	return(sigma)
	} 

robmalahanobis <- apply(((t(t(data)-center))%*%inver)*(t(t(data)-center)),1,sum) # robust malahanobis ditance
ordermalahanobis <- sort(robmalahanobis)

# comparing theoretical and empirical distribution function

bound <- qchisq(boundq,p)			# bound from which distribution function of malahanobis distances is compared
dn <- max((pchisq(ordermalahanobis,df=p)-rank(ordermalahanobis)/n+1/n)*(ordermalahanobis>=bound)) # maximum positiv! difference between theoretical and empirical distribution function, compared from borderquantile (+1/n is result of jumpdiscontinuity of empirical distribution) 

alphan <- n-floor(dn*n)		# calculation of the cutoff value
alphan <- ordermalahanobis[alphan]
weight <- (robmalahanobis < alphan)	 # cut all observations which are greater  then the observation where the difference between the distribution functions was maximal
gooddata <- data[weight==1,]

return(cor(gooddata))
}
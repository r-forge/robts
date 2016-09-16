#################
# calculates the acf using multivariate correlation estimators
# input
# timeseries: time series without NAs as vector
# lag.max: maximl lag of interest
# method: multivariate correlation estimator (weighted MCD, raw MCD, fully efficient MCD, S, Tyler-M and Stahel-Donoho) are available
# ...: further arguments for the correlation estimators
#################
 

acfmulti <- function(x, lag.max, multi.method="weightedMCD", ...) {
  timeseries <- x
  n <- length(timeseries)
  lags <- 1:lag.max
  
  # choosing the multivariate estimator:
  if (multi.method=="weightedMCD") corest <- function(x) covMcd(x)$cov
  if (multi.method=="rawMCD") corest <- function(x) covMcd(x)$raw.cov
  if (multi.method=="Stahel-Donoho") corest <- function(x) CovSde(x)@cov
  if (multi.method=="S") corest <- function(x) CovSest(x, ...)@cov
  if (multi.method=="M") corest <- function(x) mvhuberM(x, ...)$scatter
  if (multi.method=="reweight") corest <- function(x) Corefw(x, ...)
  if (multi.method=="Tyler") corest <- function(x) tyler.shape(x, location = rep(median(timeseries), lag.max + 1))
  if (multi.method=="sscor") corest <- function(x) sscor(x, location=rep(median(timeseries), lag.max+1), standardized=FALSE, pdim=TRUE)
  if (!any(multi.method==c("weightedMCD", "rawMCD", "Stahel-Donoho", "S", "reweight", "Tyler", "M","sscor"))) {
  	warning("This is no suitable correlation estimator. The weighted MCD is used instead.")
  	corest <- function(x) covMcd(x)$cov
  }
  
  corestimator <- function(x) {
  	A <- try(corest(x), silent = TRUE)
  	if (inherits(A, "try-error")) {
  		warning("Calculation of the estimator failed.")
  		return(NA)
  	}
  	return(A)
  }
  
  # building the multidimensional data matrix (insert 0 where no data is available to make it comparable to the empirical acf):
  A <- matrix(0, ncol=lag.max+1, nrow=n+lag.max)
  A[1:n,1] <- timeseries
  for (i in lags) {
  A[(i+1):(n+i),(i+1)] <- timeseries
  }
  
  # calculating the acf:
  acfvalues <- numeric(length(lags))
  covmatrix <- corestimator(A)
  if(length(covmatrix)==1) return(NA) # if something has failed before
  cormatrix <- cov2cor(covmatrix)
  for (i in lags) {
  acfvalues[i] <- mean(offdiag(cormatrix, i))  # averaging the off-diagonals
  }
  
  return(acfvalues)
}

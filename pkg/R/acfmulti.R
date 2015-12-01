#################
# calculates the acf using multivariate correlation estimators
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximl lag of interest
# method: multivariate correlation estimator (weighted MCD, raw MCD, fully efficient MCD, S, Tyler-M and Stahel-Donoho) are available
# ...: further arguments for the correlation estimators
#################
 

acfmulti <- function(timeseries,maxlag,multi.method="weightedMCD",...) {

n <- length(timeseries)


# choosing the right multivariate estimator

if (multi.method=="weightedMCD") corest <- function(x) covMcd(x)$cov
if (multi.method=="rawMCD") corest <- function(x) covMcd(x)$raw.cov
if (multi.method=="Stahel-Donoho") corest <- function(x) CovSde(x)@cov
if (multi.method=="S") corest <- function(x) CovSest(x, ...)@cov
if (multi.method=="M") corest <- function(x) mvhuberM(x, ...)$scatter
if (multi.method=="reweight") corest <- function(x) Corefw(x, ...)
if (multi.method=="Tyler") corest <- function(x) tyler.shape(x, location = rep(median(timeseries), maxlag + 1))
if (multi.method=="sscor") corest <- function(x) sscor(x,location=rep(median(timeseries),maxlag+1),standardized=FALSE,pdim=TRUE)

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



# building the multidimensional datamatrix (insert 0 where no data is available to make it comparable to the empirical acf)

A <- matrix(0,ncol=maxlag+1,nrow=n+maxlag)
A[1:n,1] <- timeseries

for (i in 1:maxlag) {
A[(i+1):(n+i),(i+1)] <- timeseries
}

# calculating the acf

acfres <- numeric(maxlag)
covmatrix <- corestimator(A)
if(length(covmatrix)==1) return(NA) # is something failed before
cormatrix <- cov2cor(covmatrix)
for (i in 1:maxlag) {
acfres[i] <- mean(offdiag(cormatrix,i))  # averaging the off-diagonals
}
return(acfres)
}

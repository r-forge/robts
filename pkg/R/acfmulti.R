#################
# calculates the acf using multivariate correlation estimators
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximl lag of interest
# method: multivariate correlation estimator (weighted MCD, raw MCD, fully efficient MCD, S, Tyler-M and Stahel-Donoho) are available
# ...: further arguments for the correlation estimators
#################
 

acfmulti <- function(timeseries,maxlag,method="weightedMCD",...) {

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

# choosing the right multivariate estimator

if (method=="weightedMCD") {
	corestimator <- function(x) {
		A <- try(covMcd(x)$cov,silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("calculation of the MCD failed")
			return(NA)
			}
		return(A)
		}
	}

if (method=="rawMCD") {
	corestimator <- function(x) {
		A <- try(covMcd(x)$raw.cov,silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("calculation of the MCD failed")			
			return(NA)
			}
		return(A)
		}
	}

if (method=="Stahel-Donoho") {
	corestimator <- function(x) {
		A <- try(CovSde(x)@cov,silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("Calculation of the Stahel-Donoho-estimator failed.")
			return(NA)
			}
		return(A)
		}
	}

if (method=="S") {
	corestimator <- function(x) {
		A <- try(CovSest(x,...)@cov,silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("Calculation of the S-estimator failed")
			return(NA)
			}
		return(A)
		}
	}

if (method=="M") {
	corestimator <- function(x) {
		A <- try(mvhuberM(x,...)$scatter,silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("Calculation of the M-estimator failed")
			return(NA)
			}
		return(A)
		}
	}

if (method=="reweight") {
	corestimator <- function(x) {
		A <- try(Corefw(x,...),silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("calculation of the fully efficient reweighted estimator failed.")
			return(NA)
			}
		return(A)
		}
	}

if (method=="Tyler") {
	Median <- median(timeseries)
	corestimator <- function(x) {
		Erg <-  try(tyler.shape(x,location=rep(Median,maxlag+1)),silent=TRUE)
		if (inherits(Erg,"try-error")){
			warning("Estimation of Tyler-shape-matrix failed.")
			return(NA)
			}
		return(Erg)
		}
	}

if (!any(method==c("weightedMCD","rawMCD","Stahel-Donoho","S","reweight","Tyler","M"))) {
	warning("This is no suitable correlation estimator. The weighted MCD is used instead.")
	corestimator <- function(x) {
		A <- try(covMcd(x)$cov,silent=TRUE)
		if (inherits(A,"try-error")) {
			warning("calculation of the MCD failed")
			return(NA)
			}
		return(A)
		}
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
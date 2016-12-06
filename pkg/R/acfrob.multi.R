#################
# calculates the acf using multivariate correlation estimators
# input
# x: time series without NAs as vector
# lag.max: maximl lag of interest
# method: multivariate correlation estimator (weighted MCD, raw MCD, fully efficient MCD, S, Tyler-M and Stahel-Donoho) are available
# ...: further arguments for the correlation estimators
#################
 

acfrob.multi <- function(x, lag.max, multi.method=c("weightedMCD", "rawMCD", "Stahel-Donoho", "S", "reweight", "Tyler", "M", "sscor"), ...) {
  multi.method <- match.arg(multi.method)
  n <- length(x)
  lags <- 1:lag.max
  
  # choosing the multivariate estimator:
  if (multi.method=="weightedMCD") covfn <- function(y) covMcd(y)$cov
  if (multi.method=="rawMCD") covfn <- function(y) covMcd(y)$raw.cov
  if (multi.method=="Stahel-Donoho") covfn <- function(y) CovSde(y)@cov
  if (multi.method=="S") covfn <- function(y) CovSest(y, ...)@cov
  if (multi.method=="M") covfn <- function(y) mvhuberM(y, ...)$scatter
  if (multi.method=="reweight") covfn <- function(y) Corefw(y, ...)
  if (multi.method=="Tyler") covfn <- function(y) tyler.shape(y, location = rep(median(x), lag.max + 1))
  if (multi.method=="sscor") covfn <- function(y) sscor(y, location=rep(median(x), lag.max+1), standardized=FALSE, pdim=TRUE)
  
  covestimator <- function(y) {
  	A <- try(covfn(y), silent = TRUE)
  	if (inherits(A, "try-error")) {
  		stop("Calculation of the estimator failed.")
  	}
  	return(A)
  }
  
  # building the multidimensional data matrix (insert 0 where no data is available to make it comparable to the empirical acf):
  A <- matrix(0, ncol=lag.max+1, nrow=n+lag.max)
  A[1:n,1] <- x
  for (i in lags) {
    A[(i+1):(n+i),(i+1)] <- x
  }
  
  # calculating the acf:
  acfvalues <- numeric(length(lags))
  covmatrix <- covestimator(A)
  if(length(covmatrix)==1) return(NA) # if something has failed before
  cormatrix <- cov2cor(covmatrix)
  for (i in lags) {
  acfvalues[i] <- mean(offdiag(cormatrix, i))  # averaging the off-diagonals
  }
  
 	res <- list(
   acfvalues = acfvalues,
   are = NA
  )
  	
  return(res)
}

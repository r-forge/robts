#####################
# filterrob.given : calculates a robustly filtered timeseries
# input:
#	 x: timeseries without NA as vector
#	 ar: coefficients of AR model fitted
#	 psifn: a robust Psi-Function
# scalefn: a function for robust scale estimation
# output: named list of 2:
#	filtered.ts: robustly filtered timeseries
#	residuals: residuals of AR fit
#####################


filterrob.given <- function(x, ar, psifn, locfn, scalefn, na.action = na.fail) {
  ists <- is.ts(x)
  x <- na.action(x)
  if (anyNA(x)) stop("NAs in 'x'")
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (ists) xtsp <- tsp(x)
  x <- as.vector(x)
  if(missing(psifn)) psifn <- function(x) M_psi(x, type="smooth")
  if(missing(locfn)) locfn <- function(x) scaleTau2(x, mu.too = TRUE)[1]
  if(missing(scalefn)) scalefn <- scaleTau2
  n <- length(x)
  p <- length(ar)
  
  # state transition matrix:
  Phi <- matrix(data = 0, nrow = p, ncol = p)
  Phi[1, ] <- ar
  if (p > 1) diag(Phi[2:p, 1:(p - 1)]) <- 1
  
  mu <- locfn(x)
  xhat <- numeric(n+p)  # robustly predicted values		
  ug <- numeric(n)		  # residuals
  sigmau <- scalefn(x)^2
  
  if(p==0) {
    residuals <- x - mu
    filtered <- psifn(residuals/sqrt(sigmau))
  } else {
    # Estimation of a start value for the recursion for filtering the error covariance:
    datamatrix <- matrix(nrow=n-p+1,ncol=p)
    for (i in 1:p) {
    	datamatrix[,i] <- x[(p-i+1):(n-i+1)] 
   	}
    P <- try(covMcd(datamatrix)$cov, silent=TRUE)
    if (inherits(P, "try-error")) {
    	warning("Start estimation of covariance filtering errors failed.")
    	return(NA)	
   	}
    
    # calculation of robustly filteres values:
    for(i in (p+1):(n+p)) {
    	xhat[i] <- mu+rev(Phi[1,])%*%(xhat[(i-p):(i-1)]-mu)	# one step ahead predictions
    	ug[i-p] <- x[i-p]-xhat[i]    # residuals
      M <- Phi%*%P%*%t(Phi)        # recursion for covariance of prediction error
      M[1,1] <- M[1,1]+sigmau
      mhat <- M[,1]
      if (mhat[1]<=0) {
      	warning("Prediction error variance is less or equal 0.")
      	return(NA)
     	}
      shat <- mhat[1]^0.5
      P <- M-1/shat^2*psifn(ug[i-p]/shat)/(ug[i-p]/shat)*mhat%*%t(mhat)	# recursion for variance of filtering error
      xhat[i:(i-p+1)] <- xhat[i:(i-p+1)]+mhat/shat*psifn(ug[i-p]/shat)	# robustly filtered values
    }
  filtered <-  xhat[-(1:p)] # initialisation by zeros
  residuals <- c(rep(NA, p), ug[-(1:p)])
  }
  if (ists) {
    attr(residuals, "tsp") <- attr(filtered, "tsp") <- xtsp
    attr(residuals, "class") <- attr(filtered, "class") <- "ts"
  } 
  res <- list(filtered = filtered, residuals = residuals)
  return(res)
}

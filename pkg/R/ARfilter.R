ARfilter <- function(timeseries, p, aicpenalty=function(p) {2*p}, psi.l=2, psi.0=3) {
  n <- length(timeseries)
  aicv <- numeric(p)                                                            
  partial <- numeric(p)		# partial autocorrelations
  phiacf <- numeric(p)
  
  timeseriesalt <- matrix(nrow=n,ncol=p)
  varold <- numeric(p)
  aralt <- matrix(ncol=p,nrow=p)
  acfalt <- matrix(ncol=p,nrow=p)
  
  mu <- median(timeseries)
  timeseries <- timeseries-mu
  # calculating first autocorrelation
  
  varest <- scaleTau2(timeseries)^2
  k <- psi.l
  l <- psi.0
  a <- (2*k^2*l^2)/(k-l)^3
  b <- -(l^3+k*l^2+4*k^2*l)/(k-l)^3
  d <- (2*l^2+2*k*l+2*k^2)/(k-l)^3
  e <- -(l+k)/(k-l)^3
  
  concor <- function(c1=4.5,c2=3) {
    b <- c2*qnorm(3/4)
    corfa <- 2*((1-b^2)*pnorm(b)-b*dnorm(b)+b^2)-1
    return(corfa)
  }
  
  segam <- function(x) return(.Call("filterinit2",c(timeseries,0),x,a,b,d,e,k,l,4.5,3,concor())[n+1])
  op <- try(optimize(segam,c(-1,1)),silent=TRUE)
  if (inherits(op,"try-error")){
  	warning("Optimization failed. Partial Autocorrelation can not be computed.")
  	return(NA)
  }
  partial[1] <- op$minimum
  helppar <- partial[1]
  aralt[1,1] <- partial[1]
  varold[1] <- op$objective
  phiacf[1] <- partial[1]		# estimated acf
  aicv[1] <- log(varold[1])+aicpenalty(2)/(n-1)
  # calculation of the other autocorrelations
  timeseriesalt[,1] <- .Call("filterinit2",c(timeseries,0),partial[1],a,b,d,e,k,l,4.5,3,concor())[-(n+1)]
  
  
  if (p > 1) for (j in 2:p) {
  	segam <- function(x) return(.Call("filter2",c(timeseries,0),x,varold[j-1],helppar,a,b,d,e,k,l,4.5,3,concor(),varest*ARMAacf(helppar))[n+1])
  	op <- try(optimize(segam,c(-1,1)),silent=TRUE)
  	if (inherits(op,"try-error")){
  		warning("Optimization failed. Partial Autocorrelation can not be computed.")
  		return(partial)
  	}	
  	partial[j] <- op$minimum
  	timeseriesalt[,j] <- try(.Call("filter2",c(timeseries,0),partial[j],varold[j-1],helppar,a,b,d,e,k,l,4.5,3,concor(),varest*phiacf[1:(j-1)])[1:n],silent=TRUE)
  	if (inherits(timeseriesalt[,j],"try-error")){
  		warning("Calculation failed, partial correlation can not be computed.")
  		return(NA)
  	}	
  	varold[j] <- op$objective
  	Phi <- numeric(j)
  	Phi[j] <- partial[j]
  	
  	# updating the memory j predictors
  	
  	for (i in 1:(j-1)) {
  		Phi[i] <- helppar[i]-partial[j]*helppar[j-i]
 		}
  	helppar <- Phi
  	aicv[j] <- log(varold[j])+aicpenalty(j+1)/(n-j)	
  	aralt[j,1:j] <- helppar 
  }
  
  erg <- list(partial, varold, ARMAacf(helppar)[-1], aicv, timeseriesalt+mu, aralt)
  names(erg) <- c("partial autocorrelations", "variance of innovations", "autocorrelation", "aic", "robustly filtered timeseries", "AR-coefficients")
  return(erg)
}

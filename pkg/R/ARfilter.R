ARfilter <- function(x, order.max, aicpenalty = function(p) {2*p}, psi.l = 2, psi.0 = 3) {
  n <- length(x)
  
  # create empty objects for the output:                                                            
  partial <- rep(NA, order.max)
  scale_vector <- rep(NA, order.max+1)
  aic_vector <- rep(NA, order.max+1) 
  filtcent_matrix <- matrix(nrow=n, ncol=order.max+1)
  acfvalues <- rep(NA, order.max)
  ar_list <- vector("list", order.max+1)
  
  # centering the time series:
  mu <- median(x)
  x <- x - mu
  
  # calculating first autocorrelation:
  varest <- scaleTau2(x)^2
  k <- psi.l
  l <- psi.0
  a <- (2*k^2*l^2)/(k-l)^3
  b <- -(l^3+k*l^2+4*k^2*l)/(k-l)^3
  d <- (2*l^2+2*k*l+2*k^2)/(k-l)^3
  e <- -(l+k)/(k-l)^3
  
  # consistency correction:
  conscorr <- function(c1=4.5, c2=3) {
    b <- c2*qnorm(3/4)
    corfa <- 2*((1-b^2)*pnorm(b)-b*dnorm(b)+b^2)-1
    return(corfa)
  }
  
  # AR(0):
  scale_vector[1] <- sqrt(varest)
  aic_vector[1] <- log(scale_vector[1]) + aicpenalty(1)/n
  filtcent_matrix[, 1] <- filterrob.given(x, ar=NULL)$filtered
  
  # AR(p) for p>0:
  if (order.max > 0) {
    segam <- function(z) .Call("filterinit2", c(x,0), z, a, b, d, e, k, l, 4.5, 3, conscorr())[n+1]
    op <- try(optimize(segam, c(-1,1)), silent=TRUE)
    if (inherits(op, "try-error")){
    	warning("Optimization failed for a model of order zero")  	
   		res <- list(pacf=partial, var=scale_vector^2, acf=rep(NA, order.max), aic=aic_vector, filtered=filtcent_matrix+mu, ar=ar_list)
   		return(res)
    }
    partial[1] <- op$minimum
    phi_temp <- partial[1]
    ar_list[[1+1]] <- phi_temp
    scale_vector[1+1] <- op$objective
    aic_vector[1+1] <- log(scale_vector[1+1]) + aicpenalty(2)/(n-1)
    # calculation of the other autocorrelations
    filtcent_matrix[, 1+1] <- .Call("filterinit2", c(x,0), partial[1], a, b, d, e, k, l, 4.5, 3, conscorr())[-(n+1)]
    acfvalues <- ARMAacf(phi_temp, lag.max=order.max)[-1]
  }
  
  if (order.max > 1) {
    for (p in 2:order.max) {
    	segam <- function(z) .Call("filter2", c(x,0), z, scale_vector[p+1-1], phi_temp, a, b, d, e,k, l, 4.5, 3, conscorr(), varest*acfvalues[1:(p-1)])[n+1]
    	op <- try(optimize(segam, c(-1,1)), silent=TRUE)
    	if (inherits(op, "try-error")){
    		warning(paste("Optimization failed for a model of order", p))
    		res <- list(pacf=partial, var=scale_vector^2, acf=rep(NA, order.max), aic=aic_vector, filtered=filtcent_matrix+mu, ar=ar_list)
    		return(res)
    	}	
    	partial[p] <- op$minimum
    	filtcent_matrix[, p+1] <- .Call("filter2", c(x,0), partial[p], scale_vector[p+1-1], phi_temp, a, b, d, e, k, l, 4.5, 3, conscorr(), varest*acfvalues[1:(p-1)])[1:n]
    	scale_vector[p+1] <- op$objective
    	aic_vector[p+1] <- log(scale_vector[p]) + aicpenalty(p+1)/(n-p)
    	Phi <- numeric(p)
    	Phi[p] <- partial[p]	
    	# updating the memory p predictors:
    	for (i in 1:(p-1)) {
    		Phi[i] <- phi_temp[i] - partial[p] * phi_temp[p-i]
   		}
    	phi_temp <- Phi
    	ar_list[[p+1]] <- phi_temp
      acfvalues <- ARMAacf(phi_temp, lag.max=order.max)[-1] 
    }
  }
  
  res <- list(pacf=partial, var=scale_vector^2, acf=acfvalues, aic=aic_vector, filtered=filtcent_matrix+mu, ar=ar_list)
  return(res)
}

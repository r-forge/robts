###################################
# calculates AR-fits using Generalized M estimates (see Robust Statistics, Maronna et al. chapter 8) of increasing order (up to order.max) and returns robust aic values
# input
# x: time series of interest as vector
# order.max: maximal order of AR-process
# maxit: maximal number of iterations for iterative weighting procedures of M-estimators
# epsilon: accuracy of iterated solution of M-estimators
# k1: tuning parameter for huber weights
# k2: tuning parameter for regressor
# output
# phi_matrix: matrix with fitted AR-parameters
# 	 AR modell of order p in row p+1, AR parameter s in column s
# RAICs:	robust aic-value of estimated AR processes (in ascending order)
#	smallest value in first element => white noise
#####################################

arrob.gm <- function(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p, na.action = na.fail, series = deparse(substitute(x)), maxit=10^3, delta=1/2, epsilon=10^(-4), k1=1.37, k2=1) {
  cl <- match.call()
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
  x <- na.action(as.ts(x))
  if (anyNA(x)) stop("NAs in 'x'")
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (ists) xtsp <- tsp(x)
  xfreq <- frequency(x)
  #x <- as.vector(x)
  n <- length(x)
  if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10)))) 
  if (order.max < 0L) stop("'order.max' must be >= 0")
  if (order.max >= n) stop("Argument 'order.max' must be lower than the length of the time series")
	if (!is.null(order.max)) if (order.max >= floor((n - 1) / 2)) {
		warning("Not enough data for chosen model order 'order.max'. The largest possible value of 'order.max' is used.")
		order.max <- floor((n - 1) / 2) - 1
	}
	if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10))))
	if (order.max < 1) stop("Model order 'order.max' must be greater than zero.") 
  # is delta valid?
  if (delta <= 0) {delta <- 0.01
    warning("Argument 'delta' <= 0 is not possible. Delta is set to 0.01.")
  }
  if (delta >= 1) {delta <- 0.5
    warning("Argument 'delta' >= 1 is not possible. Delta is set to 0.5.")
  }
  
  # calculating consistency correction
  kon <- conscorr(delta)[[1]]
  sd.pred_vector <- numeric(order.max+1)
  
  resid_matrix <- matrix(NA, ncol=order.max+1, nrow=n)
  phi_matrix <- matrix(NA, ncol=order.max, nrow=order.max)
  phiacf <- numeric(order.max-1)
  RAICs <- rep(NA, order.max+1)
 	colnames(phi_matrix) <- paste("phi", 1L:order.max, sep="_")
 	rownames(phi_matrix) <- 1L:order.max
  names(RAICs) <- colnames(resid_matrix) <- 0L:order.max
  
  # AR(0) process
  erg <- simultM_location(x, delta=delta, maxit=maxit, epsilon=epsilon, k1=k1, kon=kon)
  x.mean <- erg$location
  sd.pred_vector[1] <- erg$scale
  resid_matrix[, 1] <- erg$residuals*sd.pred_vector[1]
  if (sd.pred_vector[1]==0) {
    stop("Estimated variance is 0. Cannot fit any model")
	}
  RAICs[1] <- log(sd.pred_vector[1]^2) + aicpenalty(1)/n
  x <- x - x.mean	# centering
  
  # AR(1) process
  weightx <- M_wgt(x[-n]/sd.pred_vector[1], type="bisquare", k=k2)	# weights for Mallows-estimation (x dimension)
  erg <- simultM_regression(x[-1], x[-n], weightx=weightx, delta=delta, maxit=maxit, epsilon=epsilon, k1=k1, kon=kon)
  sd.pred_vector[2] <- erg$scale
  phiacf[1] <- erg$coef
  phi_matrix[1, 1] <- erg$coef
  resid_matrix[2:n, 2] <- erg$residuals*sd.pred_vector[2]
  RAICs[2] <- log(sd.pred_vector[2]^2) + aicpenalty(2)/(n-1)
  if (sd.pred_vector[2]>0) {
    
    if(order.max>1) {  
      # AR(p) process
      for (p in 2:order.max) {
        # using Durbin Levinson to get a one-dimensional regression
        y <- x[(p+1):n]
        z <- x[1:(n-p)]
        for (i in 1:(p-1)) {
          y <- y - phi_matrix[p-1,i] * x[(p+1-i):(n-i)]
          z <- z - phi_matrix[p-1,p-i]*x[(p+1-i):(n-i)]
        }
        # calculating the covariance matrix of the (multivariate) independent variables
        C <- diag(rep(1,p))
        for (i in 1:(p-1)) {
          C <- C + offdiag(rep(phiacf[i], p-i), i)
          C <- C + offdiag(rep(phiacf[i], p-i), -i)
        }
        # defining multivariate independent variables to calculate weights of Mallows-type
        xma <- matrix(ncol=p, nrow=n-p)
        for (i in 1:p) {
          xma[, i] <- x[i:(n-p+i-1)]
        }
        C <- C*sd.pred_vector[1]^2
        dt <- try(mahalanobis(xma, center=FALSE, cov=C)/p, silent=TRUE) # weights of independent variables
        if (inherits(dt, "try-error")) {
        	warning("Calculation of Mahalanobis distances failed. The acf might be not positiv definite. No AR models of higher order are considered.")
          break
        }
        if (sum(dt<0) > 0) {
        	warning("The acf is not positiv definit. No AR models of higher order are considered.")
        	break
       	}
        weightx <- M_wgt(sqrt(dt), type="bisquare", k=k2)
        erg <- simultM_regression(y, z, weightx, delta=delta, maxit=maxit, epsilon=epsilon, k1=k1, kon)
        sd.pred_vector[p+1] <- erg$scale
        resid_matrix[(p+1):n, p+1] <- erg$residuals*sd.pred_vector[p+1]
        RAICs[p+1] <- log(sd.pred_vector[p+1]^2) + aicpenalty(p+1)/(n-p)
        phi_matrix[p, p] <- erg$coef
        # updating AR-Parameter by Durbin-Levinson
        for (i in 1:(p-1)) {
          phi_matrix[p, i] <- phi_matrix[p-1, i] - erg$coef*phi_matrix[p-1, p-i]
        }
        # estimating acf for calculating the covariance matrix of the independent variables
        phiacf[p] <- sum(phi_matrix[p-1, 1:(p-1)] * phiacf[(p-1):1]) + erg$coef * (1 - sum(phi_matrix[p-1, 1:(p-1)] * phiacf[1:(p-1)]))
      }
    }
  } else {
    warning("Estimated variance of a AR model of order 1 is zero. An independence model is fitted.")
  }
  order_selected <- if(aic) which.min(RAICs)[[1]]-1 else order.max
  coeff <- if(order_selected==0) NULL else as.vector(phi_matrix[order_selected, 1:order_selected]) 
	resid_selected <- resid_matrix[, order_selected+1]	
	var.pred <- sd.pred_vector[order_selected+1]^2	
	partialacf <- diag(phi_matrix)
  if (ists) {
        attr(resid_selected, "tsp") <- xtsp
        attr(resid_selected, "class") <- "ts"
  }
  res <- list(
		order = order_selected,
		ar = coeff,
		var.pred = var.pred,
		x.mean = x.mean,
		x.intercept = NULL,
		aic = RAICs, #the function ar returns the difference of the AIC values with the lowest one 
		n.used = n,
		order.max = order.max,
		partialacf = array(partialacf, dim=c(length(partialacf), 1, 1)),
		resid = resid_selected,
		method = "GM",
		series = series,
		frequency = xfreq,
		call = cl,
		asy.var.coef = NULL
	)
	class(res) <- "ar"
	return(res)  
}



#######################################
# auxiliary function: estimates simultaneous mean and scale by M estimation
# input
# x: one dimensional random variables as vector
# delta: tunign parameter for M estimator, determines breakpoint, 0.5 is recommended
# epsilon: accuracy of iterated solution
# maxit: maximal number of iterations of iterative reweighting procedure
# k1: tuning parameter for huber weights
# kon: consistency-correction under normal distribution, can be calculated by conscorr
# output
# meanv: estimated mean
# sigv: estimated standard deviation
# residuals: x-meanv
########################################


simultM_location <- function(x, delta=1/2, maxit=10^3, epsilon=10^(-4), k1=1.37, kon=conscorr(delta)) {
  n <- length(x)
  
  # is delta valid?
  if (delta <= 0) {
    delta <- 0.01
    warning("Delta <= 0 is not possible. Delta is set to 0.01.")
  }
  if (delta >= 1) {
    delta <- 0.5
    warning("Delta >= 1 is not possible. Delta is set to 0.5.")
  }
  
  # calculating start values
  meanv <- mean(x)
  sigv <- mad(x)*kon
  
  # simultaneous M-estimation of location and scale using an iterative reweighting procedure
  for (i in (1:maxit)) {
  	if(sigv==0) {
      warning("estimated variance is 0")
     	return(c(meanv,sigv))
    }
  	resi <- (x-meanv)/sigv
  	we <- M_wgt(resi, type="huber", k=k1)		# Huber weights for location
  	meanvn <- sum(we*x)/sum(we)
  	resi <- resi/1.54764
  	we <- M_wgt(resi, type="bisquare", k=1)		# bisquare weights for scale
  	sigvn <- sqrt(sigv^2/n/delta*sum(resi^2*we))
  	if((abs(meanvn-meanv)<epsilon*sigv)&(abs(sigvn/sigv-1)<epsilon)) break	# stopping rule
  	sigv <- sigvn
  	meanv <- meanvn
  	if(i==maxit) warning("Iteration may not be converged")
  	}
  erg <- list(location=meanvn, scale=sigvn/kon, residuals=resi)
  return(erg)
}


#######################################
# auxiliary function: estimates simultaneous one dimensional regressioncoefficient and scale by M estimation
# input
# x: one dimensional independent variable as vector
# y: one dimensional dependent variable as vector
# weightx: weights (depending on Malahanobisdistances of x), generated by bestAR
# delta: tunign parameter for M estimator, determines breakpoint, 0.5 is recommended
# epsilon: accuracy of iterated solution
# maxit: maximal number of iterations of iterative reweighting procedure
# k1: tuning parameter for huber weights
# kon: consistency-correction under normal distribution, can be calculated by conscorr
# output
# beta: regressio coefficient
# sigv: estimated standard deviation
########################################


simultM_regression <- function(y, x, weightx, delta=1/2, maxit=10^3, epsilon=10^(-4), k1=1.37,kon=conscorr(delta)) {
  n <- length(x)
  
  # is delta valid?
  if (delta <= 0) {delta <- 0.01
  warning("Delta <= 0 is not possible. Delta is set to 0.01.")
  }
  if (delta >= 1) {delta <- 0.5
  warning("Delta >= 1 is not possible. Delta is set to 0.5.")
  }
  
  # calculate start values
  erg0 <- ltsReg(y~x-1)
  beta <- erg0$coefficients
  sigv <- erg0$scale*kon
  
  # simultaneous M-estimation of regression and scale using an iterative reweighting procedure
  for (i in (1:maxit)) {
  	if(sigv==0) {warning("estimated variance is 0")
          	return(c(beta,sigv))
  		}
  	resi <- (y-x*beta)/sigv
  	weightres <- M_wgt(resi, type="huber", k=k1)	# Huberweights for Regression
  	beta <- sum(weightres*weightx*x*y)/sum(weightres*weightx*x^2)	# Mallows-Type
  	resit <- resi/1.54764
  	we <- M_wgt(resit, type="bisquare", k=1)			# bisquareweights for scale
  	sigvn <- sqrt(sigv^2/n/delta*sum(resit^2*we))
  	resi2 <- (y-x*beta)/sigvn
  if(max(abs(resi-resi2))<epsilon) break	# stopping rule
  if(i==maxit) warning("Iteration may not be converged")
  sigv <- sigvn
  resi <- resi2
  }
  erg <- list(coef=beta, scale=sigvn/kon, residuals=resi2)
  return(erg)
}


#######################################
# auxiliary function: calculates consistency corrections for scale M-estimator based on bisquare weights with tuning parameter k
# input: tuning parameter k of bisquare weights 
# output: consistency correction factor
#######################################

conscorr <- function(x) {

  # loading consistency corrections (first column: delta, second column: consistency factor)
  corvalues <- get(load(system.file("extdata", "deltacorrection", package = "robts")))
  n <- length(corvalues[,1])
  
  # if delta is not in the usual interval
  
  if (x<corvalues[1,1]) {
  	warning("delta is to small, variance estimation will not be consistent")
  	kon <- corvalues[1,2]
  	}
  if (x==corvalues[1,1]) kon <- corvalues[1,2]
  if (x==corvalues[n,1]) kon <- corvalues[n,2]
  if (x>corvalues[n,1]) {
  	warning("delta is very large, variance estimation might not be consistent")
  	kon <- corvalues[n,2]
  	}
  
  # linear Interpolation
  if ((corvalues[1,1] < x)&( x < corvalues[n,1])) {
  	index <- max(which(x>corvalues[,1]))
  	kon <- corvalues[index,2]+(corvalues[index+1,2]-corvalues[index,2])/(corvalues[index+1,1]-corvalues[index,1])*(x-corvalues[index,1])
  	}
  return(kon)
}
 



##################
# calculating the acf using partial autokorrelations and rank estimators
# input
# x: time series without NAs as vector
# lag.max: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


acfrob.partrank <- function(x, lag.max, cor.method=c("spearman", "kendall", "quadrant", "gaussian", "masarotto"), biascorr = TRUE, partial = FALSE) {
  n <- length(x)
  cor.method <- match.arg(cor.method)

  # choosing the right estimator:  
  if (cor.method=="gaussian") {
  	correlation <- function(x, y, biascorr) {
  		n <- length(x)
  		
  		# calculating the consistency factor:  
  		i <- 1:n
  		cn <- 1/sum(qnorm(i/(n+1))^2)
  
  		# calculating the correlation:  	
  		xRank <- rank(x)
  		yRank <- rank(y)
  		result <- cn * qnorm(xRank/(n+1)) %*% qnorm(yRank/(n+1))
  		return(result)
		}
	}
  if (cor.method=="spearman")  {
  	correlation <- function(x, y, biascorr) {
 	    corvalue_biased <- cor(x, y, method = "spearman")
 	    if(!biascorr) return(corvalue_biased)
      corvalue <- 2*sin(corvalue_biased/6*pi)
      return(corvalue)
    }
 	}
  if (cor.method=="kendall")  {
   	correlation <- function(x, y, biascorr) {
 	    corvalue_biased <- cor(x, y, method = "kendall")
 	    if(!biascorr) return(corvalue_biased)
      corvalue <- sin(corvalue_biased*pi/2)
      return(corvalue)
    }
 	}
  if (cor.method=="quadrant") {
  	globalmedian <- median(x)	
  	correlation <- function(x, y, biascorr) {
  		x <- sign(x-globalmedian)
  		y <- sign(y-globalmedian)
  		n <- length(x)
  		corvalue_biased <- (t(x)%*%y)/n
  		if(!biascorr) return(corvalue_biased)
  		corvalue <- sin(corvalue_biased*pi/2)
  		return(corvalue)
 		}
 	}
  if (cor.method=="masarotto") {
    x <- x - median(x) # centering time series
  	correlation <- function(x, y, biascorr) BurgM(x, y)
  } 
  
  # calculating partial autocorrelations:  
  a <- matrix(ncol=lag.max, nrow=lag.max)	# to save changing auxiliary parameters
  phi <- numeric(lag.max)	# partial autocorrelations
  rho <- numeric(lag.max)	# autocorrelations
  
  # starting the recursion:  
  a[1,1] <- correlation(x[-1], x[-n], biascorr=biascorr)	
  phi[1] <- a[1,1]
  rho[1] <- a[1,1]
  
  # higher autocorrelations:  
  for (H in 2:lag.max) {
  	uH <- x[(H+1):n]	# foreward residuals
  	for (i in 1:(H-1)) {
  		uH <- uH - a[H-1, i] * x[(H+1-i):(n-i)]
  		}
  	vH <- x[1:(n-H)]	# backward residuals
  	for (i in 1:(H-1)) {
  		vH <- vH - a[H-1, i] * x[(1+i):(n-H+i)]
  		}
  	phi[H] <- correlation(uH, vH, biascorr=biascorr) # partial autocorrelation
  	a[H,H] <- phi[H]
  	for (i in 1:(H-1))    {
  		a[H,i] <- a[H-1, i] - phi[H] * a[H-1, (H-i)] # updating parameters
  		}  
  	rho[H] <- a[H-1, 1:(H-1)] %*% rho[(H-1):1] + phi[H] * (1-a[H-1, 1:(H-1)] %*% rho[1:(H-1)])	# calculating the autocorrelation
 	}
 	acfvalues <- if(partial) phi else rho
  
 	res <- list(
   acfvalues = acfvalues,
   are = NA
  )
  	
  return(res)
}

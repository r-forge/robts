##################
# calculating the acf using rank estimators
# input
# x: time series without NAs as vector
# lag.max: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


acrob.frank <- function(x, lag.max, cor.method=c("gaussian", "spearman", "kendall", "quadrant", "masarotto"), biascorr = TRUE) {
  n <- length(x)
  lags <- 1:lag.max
  cor.method <- match.arg(cor.method)
  
  # choosing the right estimator:  
  if (cor.method=="gaussian") {
  	# transformation into normal scores:
  	x_transformed <- qnorm(rank(x)/(n+1))
  
  	# calculating the consistency factor:
  	cn <- sum(qnorm((1:n)/(n+1))^2)
  	
  	acfvalues <- acf(x_transformed, lag.max=lag.max, plot=FALSE, type="cov")$acf[-1]/cn*n
  	return(acfvalues)
 	}
  
  if (cor.method=="spearman") {
  	# transformation into centred ranks:
  	meanrank <- (n+1)/2 # mean rank
  	x_transformed <- rank(x) - meanrank
      
  	# calculating the consistency factor:
  	cn <- sum((1:n)^2) - n*meanrank^2
   	
  	acfvalues_biased <- acf(x_transformed, lag.max=lag.max, type="cov", plot=FALSE, demean=FALSE)$acf[-1]/cn*n
  	if(!biascorr) return(acfvalues_biased)
    
    acfvalues <- 2*sin(acfvalues_biased*pi/6)
    return(acfvalues)
 	}
  
  # for the other methods a function corestimation ist defined and applied to the time series in a second step:
  
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
  
  # calculation of the acf:  
  acfvalues <- numeric(length(lags))
  for (i in lags) {
  	acfvalues[i] <- correlation(x[1:(n-i)], x[(i+1):n], biascorr=biascorr)
 	}
  return(acfvalues)
}

##################
# calculating the acf using rank estimators
# input
# x: time series without NAs as vector
# lag.max: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


acfrank <- function(x, lag.max, rank.method=c("gaussian", "spearman", "kendall", "quadrant", "masarotto")) {
  n <- length(x)
  lags <- 1:lag.max
  rank.method <- match.arg(rank.method)
  
  # centering x for masarotto approach:  
  if (rank.method=="masarotto") x <- x - median(x)
    
  # choosing the right estimator
  
    if (rank.method=="gaussian") {
  	# transformation into normal scores:
  	x_transformed <- qnorm(rank(x)/(n+1))
  
  	# calculating the consistency factor:
  	cn <- sum(qnorm((1:n)/(n+1))^2)
  	
  	acv <- acf(x_transformed, lag.max=lag.max, plot=FALSE, type="cov")$acf[-1]
  	return(acv/cn*n)
  	}
  
  if (rank.method=="spearman") {
  	# transformation into centred ranks:
  	meanrank <- (n+1)/2 # mean rank
  	x_transformed <- rank(x) - meanrank
      
  	# calculating the consistency factor:
  	cn <- sum((1:n)^2) - n*meanrank^2
   	
  	acv <- acf(x_transformed, lag.max=lag.max, type="cov", plot=FALSE, demean=FALSE)$acf[-1]
  	return(2*sin(acv/cn*n*pi/6))
  	}
  
  # for the other methods a function corestimation ist defined and applied to the time series in a second step:
  
  if (rank.method=="kendall") {
  	corestimation <- function(x, y) {sin(cor(x, y, method="kendall")*pi/2)}
  	}
  
  if (rank.method=="quadrant") {
  	globalmedian <- median(x)	
  	corestimation <- function(x, y) {
  		x <- sign(x-globalmedian)
  		y <- sign(y-globalmedian)
  		n <- length(x)
  		erg <- (t(x)%*%y)/n
  		return(sin(erg*pi/2))
  		}
  	}
  if (rank.method=="masarotto") {
  	corestimation <- function(x, y) BurgM(x, y)
  }
  
  # calculation of the acf:  
  acfvalues <- numeric(length(lags))
  for (i in lags) {
  	acfvalues[i] <- corestimation(x[1:(n-i)], x[(i+1):n])
 	}
  return(acfvalues)
}

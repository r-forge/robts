##################
# calculating the acf using partial autokorrelations and rank estimators
# input
# x: time series without NAs as vector
# lag.max: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


acfpartrank <- function(x, lag.max, partial.method=c("gaussian", "spearman", "kendall", "quadrant", "masarotto")) {
  partial.method <- match.arg(partial.method)
  timeseries <- x
  n <- length(timeseries)
 
  # centering timeseries for Masarotto approach:
  if (partial.method=="masarotto") timeseries <- timeseries - median(timeseries)
  
  # choosing the right estimator:  
  if (partial.method=="gaussian") {
  	correlation <- function(x, y) {
  		n <- length(x)
  		
  		# calculating the consistency factor:  
  		i <- 1:n
  		cn <- 1/sum(qnorm(i/(n+1))^2)
  
  		# calculating the correlation:  	
  		xRang <- rank(x)
  		yRang <- rank(y)
  		Kor <- cn * qnorm(xRang/(n+1)) %*% qnorm(yRang/(n+1))
  		return(Kor)
		}
	}
  if (partial.method=="spearman")  {
  	correlation <- function(x, y) {2*sin(cor(x,y,method="spearman")/6*pi)}
  	}
  if (partial.method=="kendall")  {
  	correlation <- function(x, y) {sin(cor(x,y,method="kendall")*pi/2)}
  	}
  if (partial.method=="quadrant") {
  	Median <- median(timeseries)	
  	correlation <- function(x, y) {
  		x <- sign(x-Median)
  		y <- sign(y-Median)
  		n <- length(x)
  		erg <- (t(x)%*%y)/n
  		return(sin(erg*pi/2))
  		}
  	}
  if (partial.method=="masarotto") {
  	correlation <- function(x,y) BurgM(x, y)
  } 
  
  # calculating partial autocorrelations:  
  a <- matrix(ncol=lag.max, nrow=lag.max)	# to save changing auxiliary parameters
  phi <- numeric(lag.max)	# partial autocorrelations
  rho <- numeric(lag.max)	# autocorrelations
  
  # starting the recursion:  
  a[1,1] <- correlation(timeseries[-1], timeseries[-n])	
  phi[1] <- a[1,1]
  rho[1] <- a[1,1]
  
  # higher autocorrelations:  
  for (H in 2:lag.max) {
  	uH <- timeseries[(H+1):n]	# foreward residuals
  	for (i in 1:(H-1)) {
  		uH <- uH - a[H-1, i] * timeseries[(H+1-i):(n-i)]
  		}
  	vH <- timeseries[1:(n-H)]	# backward residuals
  	for (i in 1:(H-1)) {
  		vH <- vH - a[H-1, i] * timeseries[(1+i):(n-H+i)]
  		}
  	phi[H] <- correlation(uH, vH) # partial autocorrelation
  	a[H,H] <- phi[H]
  	for (i in 1:(H-1))    {
  		a[H,i] <- a[H-1, i] - phi[H] * a[H-1, (H-i)] # updating parameters
  		}  
  	rho[H] <- a[H-1, 1:(H-1)] %*% rho[(H-1):1] + phi[H] * (1-a[H-1, 1:(H-1)] %*% rho[1:(H-1)])	# calculating the autocorrelation
 	}

  return(rho)
}

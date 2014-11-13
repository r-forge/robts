##################
# calculating the acf using rank estimators
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


rankacf <- function(timeseries,maxlag,method="gaussian") {

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

# centering timeseries for masarotto approach

if (method=="masarotto") timeseries <- timeseries - median(timeseries)


# choosing the right estimator


if (method=="gaussian") {
	# transformation into ranks
	rang <- rank(timeseries)
	n <- length(timeseries)

	# calculating the consistency factor
	i <- 1:n
	cn <- sum(qnorm(i/(n+1))^2)

	kor <- numeric(maxlag)
	for (i in 1:maxlag) {
		kor[i] <- sum(qnorm(rang[1:(n-i)]/(n+1))*qnorm(rang[(i+1):n]/(n+1)))/cn
		}
	return(kor)
	}

if (method=="spearman") {
	# transformation into ranks
	rang <- rank(timeseries)
	n <- length(timeseries)

	# calculating mean value
	mv <- (n+1)/2

	# calculating the consistency factor
	i <- 1:n
	cn <- sum(i^2)-n*mv^2
	
	
	

	kor <- numeric(maxlag)
	for (i in 1:maxlag) {
		kor[i] <- sum((rang[1:(n-i)]-mv)*(rang[(i+1):n]-mv))/cn
		}
	kor <- 2*sin(kor*pi/6)
	return(kor)
	}

if (method=="kendall") {
	corestimation <- function(x,y) {sin(cor(x,y,method="kendall")*pi/2)}
	}

if (method=="quadrant") {
	Median <- median(timeseries)	
	corestimation <- function(x,y) {
		x <- sign(x-Median)
		y <- sign(y-Median)
		n <- length(x)
		erg <- (t(x)%*%y)/n
		return(sin(erg*pi/2))
		}
	}
if (method=="masarotto") {
	corestimation <- function(x,y) masarotto(x,y)
}


if (!any(method==c("gaussian","spearman","kendall","quadrant","masarotto"))) {
	warning("This is no suitable correlation estimator. Kendalls-Tau is used instead.")
	corestimation <- function(x,y) {sin(cor(x,y,method="kendall")*pi/2)}
	}

# calculation of the acf

acfvalues <- numeric(maxlag)
for (i in 1:maxlag) {
	acfvalues[i] <- corestimation(timeseries[1:(n-i)],timeseries[(i+1):n])
	}
return(acfvalues)
}
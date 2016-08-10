##################
# calculating the acf using rank estimators
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


acfrank <- function(timeseries,maxlag,rank.method="gaussian") {

n <- length(timeseries)

# centering timeseries for masarotto approach

if (rank.method=="masarotto") timeseries <- timeseries - median(timeseries)


# choosing the right estimator


if (rank.method=="gaussian") {
	# transformation into ranks
	rang <- rank(timeseries)
	n <- length(timeseries)

	# calculating the consistency factor
	i <- 1:n
	cn <- sum(qnorm(i/(n+1))^2)
	beob <- qnorm(rang/(n+1))
	acv <- acf(beob,lag.max=maxlag,plot=FALSE,type="cov")$acf[-1]
	return(acv/cn*n)
	}

if (rank.method=="spearman") {
	# transformation into ranks
	rang <- rank(timeseries)
	n <- length(timeseries)

	# calculating mean value
	mv <- (n+1)/2

	# calculating the consistency factor
	i <- 1:n
	cn <- sum(i^2)-n*mv^2
	
	
	acv <- acf(rang-mv,lag.max=maxlag,type="cov",plot=FALSE,demean=FALSE)$acf[-1]
	acv <- 2*sin(acv/cn*n*pi/6)
	return(acv)
	}

if (rank.method=="kendall") {
	corestimation <- function(x,y) {sin(cor(x,y,method="kendall")*pi/2)}
	}

if (rank.method=="quadrant") {
	Median <- median(timeseries)	
	corestimation <- function(x,y) {
		x <- sign(x-Median)
		y <- sign(y-Median)
		n <- length(x)
		erg <- (t(x)%*%y)/n
		return(sin(erg*pi/2))
		}
	}
if (rank.method=="masarotto") {
	corestimation <- function(x,y) BurgM(x,y)
}


if (!any(rank.method==c("gaussian","spearman","kendall","quadrant","masarotto"))) {
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

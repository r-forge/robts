##################
# calculating the acf using partial autokorrelations and rank estimators
# input
# timeseries: timeseries without NAs as vector
# maxlag: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: autocorrelation function as vector
##################


Autokor <- function(timeseries,maxlag,method="spearman") {

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
	correlation <- function(x,y) {
		n <- length(x)
		
		# calculating the consistency factor

		i <- 1:n
		cn <- 1/sum(qnorm(i/(n+1))^2)

		# calculating the correlation
	
		xRang <- rank(x)
		yRang <- rank(y)
		Kor <- cn*qnorm(xRang/(n+1))%*%qnorm(yRang/(n+1))
		return(Kor)
		}
	}

if (method=="spearman")  {
	correlation <- function(x,y) {2*sin(cor(x,y,method="spearman")/6*pi)}
	}

if (method=="kendall")  {
	correlation <- function(x,y) {sin(cor(x,y,method="kendall")*pi/2)}
	}

if (method=="quadrant") {
	Median <- median(timeseries)	
	correlation <- function(x,y) {
		x <- sign(x-Median)
		y <- sign(y-Median)
		n <- length(x)
		erg <- (t(x)%*%y)/n
		return(sin(erg*pi/2))
		}
	}
if (method=="masarotto") {
	correlation <- function(x,y) masarotto(x,y)
}

if (!any(method==c("gaussian","spearman","kendall","quadrant","masarotto"))) {
warning("hTis is no suitable correlation estimator. Gaussian-rank-correlation is used instead.")
	correlation <- function(x,y) {
		n <- length(x)
		
		# calculating the consistency factor

		i <- 1:n
		cn <- 1/sum(qnorm(i/(n+1))^2)

		# calculating the correlation
	
		xRang <- rank(x)
		yRang <- rank(y)
		Kor <- cn*qnorm(xRang/(n+1))%*%qnorm(yRang/(n+1))
		return(Kor)
		}
	}

# calculating partial autocorrelations

a <- matrix(ncol=maxlag,nrow=maxlag)		# to save changing auxiliary parameters
phi <- numeric(maxlag)				# partial autocorrelations
rho <- numeric(maxlag)				# autocorrelations

# starting the recursion

a[1,1] <- correlation(timeseries[-1],timeseries[-n])	
phi[1] <- a[1,1]
rho[1] <- a[1,1]

# higher autocorrelations

for (H in 2:maxlag) {
	uH <- timeseries[(H+1):n]	# foreward residuals
	for (i in 1:(H-1)) {
		uH <- uH-a[H-1,i]*timeseries[(H+1-i):(n-i)]
		}
	vH <- timeseries[1:(n-H)]	# backward residuals
	for (i in 1:(H-1)) {
		vH <- vH-a[H-1,i]*timeseries[(1+i):(n-H+i)]
		}
	phi[H] <- correlation(uH,vH)	# partial autocorrelation
	a[H,H] <- phi[H]

	for (i in 1:(H-1))    {
		a[H,i] <- a[H-1,i]-phi[H]*a[H-1,(H-i)]		# updating parameters
		}

	rho[H] <- a[H-1,1:(H-1)]%*%rho[(H-1):1]+phi[H]*(1-a[H-1,1:(H-1)]%*%rho[1:(H-1)])	# calculating the autocorrelation
	}

return(rho)
}
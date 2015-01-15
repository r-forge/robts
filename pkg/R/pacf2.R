##################
# calculating the partial autokorrelations (by rank estimation)
# input
# x: x without NAs as vector
# lag.max: maximal lag of interest
# method: which correlation estimator to use (Gaussian rank correlation, Spearman, Kendall, Quadrant-correlation are available)
# output: partial autocorrelation function as vector
##################


pacf2 <- function(x, lag.max = NULL, method = c("spearman", "gaussian", "kendall", "quadrant", "masarotto"),
     na.action = na.fail) {

method <- match.arg(method)

# protective measures

x <- na.action(x)

if(sum(is.na(x)) > 0) {
	warning("There are NA in your x you should use a procedure to replace this values first.")
	return(NA)
	}
n <- length(x)

if (n < 2 + lag.max) {
	warning("You can't calculate this amount of lags. We will compute the maximal possible lag of n-2")
	lag.max <- n - 2
	}
if (n < 4 * lag.max) {
	warning("Brockwell and Davis suggest to calculate only lags less n/4. Nevertheless we will calculate all lags you want, but you should be aware that the estimated acf for higher lags could be unreasonable.")
	}

# centering x for masarotto approach
if (method=="masarotto") x <- x - median(x)

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
	Median <- median(x)	
	correlation <- function(x,y) {
		x <- sign(x-Median)
		y <- sign(y-Median)
		n <- length(x)
		erg <- (t(x)%*%y)/n
		return(sin(erg*pi/2))
		}
	}
if (method=="masarotto") {
	correlation <- function(x,y) BurgM(x,y)
}

if (!any(method==c("gaussian","spearman","kendall","quadrant","masarotto"))) {
warning("This is no suitable correlation estimator. Gaussian-rank-correlation is used instead.")
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

a <- matrix(ncol=lag.max,nrow=lag.max)		# to save changing auxiliary parameters
phi <- numeric(lag.max)				# partial autocorrelations

# starting the recursion

a[1,1] <- correlation(x[-1],x[-n])	
phi[1] <- a[1,1]

# higher autocorrelations

for (H in 2:lag.max) {
	uH <- x[(H+1):n]	# foreward residuals
	for (i in 1:(H-1)) {
		uH <- uH-a[H-1,i]*x[(H+1-i):(n-i)]
		}
	vH <- x[1:(n-H)]	# backward residuals
	for (i in 1:(H-1)) {
		vH <- vH-a[H-1,i]*x[(1+i):(n-H+i)]
		}
	phi[H] <- correlation(uH,vH)	# partial autocorrelation
	a[H,H] <- phi[H]

	for (i in 1:(H-1))    {
		a[H,i] <- a[H-1,i]-phi[H]*a[H-1,(H-i)]		# updating parameters
		}
	}


	return(phi)

}
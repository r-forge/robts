## solveYW - calculation of AR model coefficients by solving the Yule-Walker equations
## input:
## 	acf: autocorrelation function
## 	p: order of the AR model
## output: AR model coefficients

solveYW <- function(acf, p) {
	p <- as.integer(p)
	stopifnot(is.numeric(acf), p >= 1)
	if (length(acf) < p) stop("lags up to order of AR model required")
	# autocorrelation function:
	acorf <- acf[1:p]
	CC <- diag(nrow = p)
	if (p > 1) for (i in 1:(p-1)) {
		CC[i, (i+1):p] <- acorf[1:(p-i)]
		CC[p-i+1, 1:(p-i)] <- acorf[(p-i):1]
	}
	ph <- solve(a = CC, b = acorf)
	return(ph)
}
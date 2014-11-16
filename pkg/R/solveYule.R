## solveYule - calculation of AR model coefficients by solving the Yule-Walker equations
## input:
## 	ts: time series
## 	acf: autocorrelation function calculated from ts
## 	p: order of the AR model
## 	robvar: should the variance be calculated robustly?
## output: AR model coefficients

solveYule <- function(ts, acf, p, robvar = TRUE) {
	p <- as.integer(p)
	stopifnot(is.numeric(acf), is.numeric(ts), p >= 1, is.logical(robvar))
	if (length(acf) < p) stop("lags up to order of AR model required")
	if (robvar) s <- scaleTau2(ts) else s <- sd(ts)
	# autocovariance function:
	acovf <- s^2 * acf[1:p]
	CC <- diag(x = s^2, nrow = p)
	if (p > 1) for (i in 1:(p-1)) {
		CC[i, (i+1):p] <- acovf[1:(p-i)]
		CC[p-i+1, 1:(p-i)] <- acovf[(p-i):1]
	}
	ph <- solve(a = CC, b = acovf)
	return(ph)
}
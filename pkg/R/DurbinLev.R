## DurbinLev - partial autocorrelation function by Durbin-Levinson algorithm
## input:
## 		ts: time series
## 		acf: autocorrelation function calculated from ts
## 		robvar: should the variance be calculated robustly?
## output: list of 2:
##		phis: list of coefficient vectors
##		nus: vector of mean squared errors

DurbinLev <- function(ts, acf, robvar = TRUE) {
	stopifnot(is.numeric(acf), is.numeric(ts), is.logical(robvar))
	if (robvar) s <- scaleTau2(ts) else s <- sd(ts)
	# autocovariance function from 0:
	gammas <- s^2 * c(1, acf)
	n <- length(gammas) - 1
	if (n < 1) stop("autocorrelation of lag 1 required at least")
	nus <- numeric(n + 1)
	phis <- vector("list", n)
	phis[[1]] <- gammas[2] / gammas[1]
	nus[1] <- gammas[1]
	nus[2] <- nus[1] * (1 - phis[[1]][1]^2)
	if(n >= 1) {
		for(i in 2:n) {
			phi_ii <- (gammas[i + 1] - sum(phis[[i - 1]] * gammas[i - 0:(i - 2)])) * 1/nus[i]
			phis[[i]] <- c(phis[[i - 1]] - phi_ii * rev(phis[[i - 1]]), phi_ii)
			nus[i + 1] <- nus[i] * (1 - phi_ii^2)
		}
	}
	res <- list(phis = phis, nus = nus)
	return(res)
}
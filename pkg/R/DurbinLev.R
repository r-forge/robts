## DurbinLev - partial autocorrelation function by Durbin-Levinson algorithm
## input:
## 		acfvalues: autocorrelation function
## output: list of 2:
##		phis: list of coefficient vectors
##		nus: vector of mean squared errors

DurbinLev <- function(acfvalues) {
	stopifnot(is.numeric(acfvalues))
	# autocorrelation function from 0:
	gammas <- c(1, acfvalues)
	n <- length(gammas) - 1
	if (n < 1) stop("autocorrelation of lag 1 required at least")
	nus <- numeric(n + 1)
	phis <- vector("list", n)
	phis[[1]] <- gammas[2]
	nus[1] <- 1
	nus[2] <- nus[1] * (1 - phis[[1]][1]^2)
	if(n >= 1) {
		for(i in 2:n) {
			phi_ii <- (gammas[i + 1] - sum(phis[[i - 1]] * gammas[i - 0:(i - 2)])) / nus[i]
			phis[[i]] <- c(phis[[i - 1]] - phi_ii * rev(phis[[i - 1]]), phi_ii)
			nus[i + 1] <- nus[i] * (1 - phi_ii^2)
		}
	}
	res <- list(phis = phis, nus = nus)
	return(res)
}
## ARopt.acf - optimally fitted AR model calculated from the acf
## Input:
##		ts: time series
##		pmax: maximum AR order considered
##		acf.fun: function for acf calculation
## Output: coefficients vector of optimal AR model


ARopt.acf <- function(tss, pmax = NULL, acf.fun = "acfmedian") {
	posfuns <- c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "trimmedcor")
	if (!any(acf.fun == posfuns)) stop("No valid ACF function.")
	tmax <- length(tss)
	lmax <- floor(tmax / 4)
	if (!is.null(pmax)) if (pmax > lmax) {
		warning("Too less acf values calculated for chosen pmax, corrected to greatest possible p.")
		pmax <- lmax
	}
	if (is.null(pmax)) pmax <- floor(min((tmax - 1) / 4, 10 * log(tmax, base = 10)))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	acf.fun <- get(acf.fun)
	acorf <- acf.fun(timeseries = tss, maxlag = lmax)
	fits <- DurbinLev(acf = acorf)
	RAICs <- numeric(pmax)
	for (p in 1:pmax) {
		D <- matrix(nrow = tmax - p, ncol = p)
		for (j in 1:p) D[, j] <- tss[(p + 1 - j):(tmax - j)]
		D <- cbind(1, D)
		ph <- fits$phis[[p]]
		ph <- c(median(tss) * (1 - sum(ph)), ph)
		resi <- tss[(p + 1):tmax] - D %*% ph
		RAICs[p] <- log(Qn(resi)^2) + 2 * p / (tmax - p)
	}
	popt <- which.min(RAICs)
	phopt <- fits$phis[[popt]]
	for (i in 1:popt) names(phopt)[i] <- paste("phi", i)
	return(phopt)
}
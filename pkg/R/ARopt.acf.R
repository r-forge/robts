## ARopt.acf - optimally fitted AR model calculated from the acf
## Input:
##		ts: time series
##		pmax: maximum AR order considered
##		acf.fun: function for acf calculation
## Output: coefficients vector of optimal AR model


ARopt.acf <- function(tss, aic = TRUE, pmax = NULL, acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim")) {
	acf.fun <- match.arg(acf.fun)
	posfuns <- c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim")
	if (!any(acf.fun == posfuns)) stop("No valid ACF function.")
	tmax <- length(tss)
	lmax <- floor(tmax / 4)
	if (!is.null(pmax)) if (pmax > lmax) {
		warning("Too less acf values calculated for chosen pmax, corrected to greatest possible p.")
		pmax <- lmax
	}
	if (is.null(pmax)) pmax <- floor(min((tmax - 1) / 4, 10 * log(tmax, base = 10)))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	acorf <- as.numeric(acfrob(tss, lag.max = lmax, fun = acf.fun, plot = FALSE)$acf)
	fits <- DurbinLev(acf = acorf)
	if (aic) {
		RAICs <- numeric(pmax + 1)
		# null model:
		resi <- tss - median(tss)
		RAICs[1] <- log(Qn(resi)^2)
		phopt <- NULL
		for (p in 1:pmax) {
			D <- matrix(nrow = tmax - p, ncol = p)
			for (j in 1:p) D[, j] <- tss[(p + 1 - j):(tmax - j)]
			D <- cbind(1, D)
			ph <- fits$phis[[p]]
			ph <- c(median(tss) * (1 - sum(ph)), ph)
			resi <- tss[(p + 1):tmax] - D %*% ph
			RAICs[p + 1] <- log(Qn(resi)^2) + 2 * p / (tmax - p)
		}
		popt <- which.min(RAICs)[1] - 1
		aic <- RAICs[popt + 1]
	} else {
		popt <- pmax
		aic <- NULL
	}
	if (popt > 0) phopt <- fits$phis[[popt]]
	return(list(coefficients = phopt, aic = aic))
}
## ARopt.YW - optimally fitted AR model calculated from the acf
## Input:
##		ts: time series
##		pmax: maximum AR order considered
##		acf.fun: function for acf calculation
## Output: coefficients of optimal AR model


ARopt.YW <- function(tss, pmax = NULL, acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"),aicpenalty=function(n,p) {return(2*p/n)}) {
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
	# null model:
	x.mean <- median(tss)
	residuals <- tss - x.mean
	var.pred <- Qn(residuals)^2
	RAICopt <- log(var.pred)
	pacfbest <- rep(0,pmax)
	phopt <- NULL
	for (p in 1:pmax) {
		D <- matrix(nrow = tmax - p, ncol = p)
		for (j in 1:p) D[, j] <- tss[(p + 1 - j):(tmax - j)]
		D <- cbind(1, D)
		ph <- solveYW(acorf, p)
		ph1 <- c(x.mean * (1 - sum(ph)), ph)
		resi <- tss[(p + 1):tmax] - D %*% ph1
		var.new <- Qn(resi)^2
		RAICnew <- log(var.new) + aicpenalty(tmax-p,p)
		if (RAICnew < RAICopt) {
			RAICopt <- RAICnew
			residuals <- resi
			var.pred <- var.new
			phopt <- ph
			pacfbest <- ARMAacf(ar=phopt,lag.max=pmax,pacf=TRUE)
		}
	}
	return(list(coefficients = phopt, aic = RAICopt,residuals=residuals,x.mean=x.mean,var.pred=var.pred,partialacf=pacfbest))
}

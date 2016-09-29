## ARopt.acf - optimally fitted AR model calculated from the acf
## Input:
##		ts: time series
##		pmax: maximum AR order considered
##		acf.approach: function for acf calculation
## Output: coefficients vector of optimal AR model 
##
##


ARopt.acf <- function(tss, pmax = NULL, acf.approach = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim","acfrobfil","acfregression"),aicpenalty=function(p) {return(2*p)},...) {
	acf.approach <- match.arg(acf.approach)
	tmax <- length(tss)
	lmax <- floor(tmax / 4)
	if (!is.null(pmax)) if (pmax > lmax) {
		warning("Too less acf values calculated for chosen pmax, corrected to greatest possible p.")
		pmax <- lmax
	}
	if (is.null(pmax)) pmax <- floor(min(c((tmax - 1) / 4, 10 * log(tmax, base = 10))))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	acorf <- as.numeric(acfrob(tss, lag.max = pmax, approach = acf.approach, plot = FALSE,...)$acf)
	if (sum(is.na(acorf))>0) {
	return(NA)
	}
	fits <- DurbinLev(acf = acorf)

		RAICs <- numeric(pmax + 1)
		# null model:
		x.mean <- median(tss)
		residuals <- tss - x.mean
		var.pred <- Qn(residuals)^2
		RAICs[1] <- log(var.pred)+aicpenalty(1)/tmax
		RAICopt <- RAICs[1]
		pacfbest <- rep(0,pmax)
		phopt <- NULL
		for (p in 1:pmax) {
			D <- matrix(nrow = tmax - p, ncol = p)
			for (j in 1:p) D[, j] <- tss[(p + 1 - j):(tmax - j)]
			D <- cbind(1, D)
			ph <- fits$phis[[p]]
			ph <- c(x.mean * (1 - sum(ph)), ph)
			resi <- tss[(p + 1):tmax] - D %*% ph
			var.new <- Qn(resi)^2
			RAICs[p+1] <- log(var.new) + aicpenalty(p+1)/(tmax-p)
			if (RAICs[p+1] < RAICopt) {RAICopt <- RAICs[p+1]
						residuals <- resi
						var.pred <- var.new
						phopt <- ph
						pacfbest <- ARMAacf(ar=phopt,lag.max=pmax,pacf=TRUE)
						}
		}
	return(list(coefficients = phopt[-1], aic = RAICs,residuals=residuals,x.mean=x.mean,var.pred=var.pred,partialacf=pacfbest))
}

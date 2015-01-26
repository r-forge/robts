## lmrobARopt - optimally fitted AR model with respect to AIC
## Input:
##		ts: time series
##		interc: with intercept?
##		##		arguments to be passed to lmrob
## Output: list of 2:
## 		coefficients: vector of coefficients (with intercept in the last place)
##		model: object of class lmrob
##		p.optimal


lmrobARopt <- function(ts, interc = TRUE, singular.ok = FALSE, pmax = NULL, ...) {
	tmax <- length(ts)
	o1 <- as.numeric(interc)
	if (!is.null(pmax)) if (pmax >= floor((tmax - 1 - o1) / 2)) {
		warning("Too less data for chosen pmax, corrected to greatest possible value.")
		pmax <- floor((tmax - 1 - o1) / 2) - 1
	}
	if (is.null(pmax)) pmax <- floor(min((tmax - 1 - o1) / 4, 10 * log(tmax, base = 10)))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	# null model:
	resi <- ts - median(ts)
	paicbest <- log(Qn(resi)^2)
	phopt <- NULL
	popt <- 0
	for (p in 1:pmax) {
		fit <- suppressWarnings(lmrobAR(ts = ts, p = p, interc = interc,
			singular.ok = singular.ok, ...)$model)
		if (fit$converged) {
			paicnew <- log(Qn(fit$residuals)^2) + 2 * p / (tmax - p)
			if (paicnew < paicbest) {
				paicbest <- paicnew
				fitbest <- fit
				popt <- p
			}
		}
	}
	if (popt == 0) {
		be <- NULL
		fitbest <- NULL
	} else {
		be <- fitbest$coefficients
		for (i in 1:popt) names(be)[i] <- paste("phi", i)
		if (interc) names(be)[popt + 1] <- "interc"
	}
	return(list(coefficients = be, model = fitbest, p.optimal = popt, aic = paicbest))
}
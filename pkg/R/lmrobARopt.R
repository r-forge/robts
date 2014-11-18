## lmrobARopt - optimally fitted AR model with respect to adjusted R-squared
## Input:
##		ts: time series
##		##		arguments to be passed to lmrob
## Output: list of 2:
## 		coefficients: vector of coefficients with intercept in the last place
##		model: object of class lmrob
##		p.optimal


lmrobARopt <- function(ts, method = "MM", singular.ok = FALSE, init = NULL) {
	tmax <- length(ts)
	if (tmax < 7) stop("Too less data for reasonable model comparison. Try p = 1.")
	pmax <- floor((tmax - 1) / 2) - 1 # AIC implausible small for utmost p
	aicbest <- +Inf
	popt <- 0
	for (p in 1:pmax) {
		fit <- suppressWarnings(lmrobAR(ts = ts, p = p, method = method, singular.ok = singular.ok, init = init)$model)
		if (fit$converged) {
			aicnew <- log(sum((fit$residuals)^2) / (tmax - p)) + 2 * p / (tmax - p)
			if (aicnew < aicbest) {
				aicbest <- aicnew
				fitbest <- fit
				popt <- p
			}
		}
	}
	if (popt == 0) stop("No convergence for any p.")
	be <- fitbest$coefficients
	for (i in 1:popt) names(be)[i] <- paste("phi", i)
	names(be)[popt + 1] <- "interc"
	return(list(coefficients = be, model = fitbest, p.optimal = popt))
}
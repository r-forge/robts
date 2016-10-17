## lmrobARopt - optimally fitted AR model with respect to AIC
## Input:
##		ts: time series
##		interc: with intercept?
##		##		arguments to be passed to lmrob
## Output: list of 2:
## 		coefficients: vector of coefficients (with intercept in the last place)
##		model: object of class lmrob
##		p.optimal


lmrobARopt <- function(ts, interc = TRUE, singular.ok = FALSE, pmax = NULL,aicpenalty=function(p) {return(2*p)}, ...) {
	tmax <- length(ts)
	o1 <- as.numeric(interc)
	if (!is.null(pmax)) if (pmax >= floor((tmax - 1 - o1) / 2)) {
		warning("Too less data for chosen pmax, corrected to greatest possible value.")
		pmax <- floor((tmax - 1 - o1) / 2) - 1
	}
	if (is.null(pmax)) pmax <- floor(min(c((tmax - 1 - o1) / 4, 10 * log(tmax, base = 10))))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	RAIC <- numeric(pmax+1)
	# null model:
	if (interc==TRUE) {	fitbest <- lmrob(ts~1)
				p.opt <- 0
				phopt <- NULL
				RAIC[1] <- log(fitbest$scale^2)+aicpenalty(1)/tmax
				RAICbest <- RAIC[1]
				resid <- fitbest$residuals
				pacfbest <- rep(0,pmax)
				var.pred <- fitbest$scale^2
				x.mean <- fitbest$coefficients
			  }	else {	x.mean <- 0
					var.pred <- (Qn(ts)^2)
					p.opt <- 0
					phopt <- NULL
					RAIC[1] <- log(var.pred)
					RAICbest <- RAIC[1]
					resid <- ts
					pacfbest <- rep(0,pmax)
				     }
			
	for (p in 1:pmax) {
		fit <- suppressWarnings(ARlmrob(ts = ts, p = p, interc = interc,
			singular.ok = singular.ok, ...)$model)
		if (fit$converged) {
			RAIC[p+1] <- log(fit$scale^2) + aicpenalty(p+interc)/(tmax-p)
			if (RAIC[p+1] < RAICbest) {
				RAICbest <- RAIC[p+1]
				fitbest <- fit
				p.opt <- p
				if (interc==TRUE) {	phopt <- fitbest$coefficients[2:(p+1)]
							x.mean <- fitbest$coefficients[1]/(1-sum(phopt))
						  }	else{	phopt <- fitbest$coefficients[1:p]}
				var.pred <- fitbest$scale^2
				resid <- as.numeric(fitbest$residuals)
				pacfbest <- ARMAacf(ar=phopt,lag.max=pmax,pacf=TRUE)
			}
		}
	}
	
	return(list(coefficients = phopt, model = fitbest, p.optimal = p.opt, aic = RAIC,x.mean=x.mean,var.pred=var.pred,partialacf=pacfbest,resid=resid))
}

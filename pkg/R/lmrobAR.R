## lmrobAR - autoregression coefficients estimated robustly via linear model
## Input:
##		ts: time series
## 		p: order of the AR model
##		interc: with intercept?
##		arguments to be passed to lmrob
## Output: list of 2:
## 		coefficients: vector of coefficients (with intercept in the last place)
##		model: object of class lmrob


lmrobAR <- function(ts, p, interc = TRUE, singular.ok = FALSE, ...) {
	p <- as.integer(p)
	stopifnot(is.numeric(ts))
	tmax <- length(ts)
	n <- tmax - p
	if (n < p + as.numeric(interc) + 1) stop("Not enough data: lesser order of AR model required.")
	y <- ts[(p + 1):tmax]
	G <- matrix(nrow = n, ncol = p)
	for (i in 1:p) G[, p + 1 - i] <- ts[i:(i + n - 1)]
	if (interc==TRUE) {	lmo <- lmrob(formula = y ~ G, singular.ok = singular.ok, ...)
				be <- c(lmo$coefficients[-1],lmo$coefficients[1])
			  }	else{	lmo <- lmrob(formula = y ~ 0 + G, singular.ok = singular.ok, ...)
					be <- lmo$coefficients
				    }		
	
	for (i in 1:p) names(be)[i] <- paste("phi", i)
	if (interc) names(be)[p + 1] <- "interc"
	return(list(coefficients = be, model = lmo))
}

## lmrobAR - autoregression coefficients estimated robustly via linear model
## Input:
##		ts: time series
## 		p: order of the AR model
##		
## Output: list of 2:
## 		coefficients: vector of coefficients with intercept in the last place
##		model: object of class lmrob
##		arguments to be passed to lmrob

lmrobAR <- function(ts, p, method = "MM", singular.ok = FALSE, init = NULL) {
	p <- as.integer(p)
	stopifnot(is.numeric(ts))
	tmax <- length(ts)
	n <- tmax - p
	if (n < p + 1) stop("Not enough data: lesser order of AR model required.")
	y <- ts[(p + 1):tmax]
	G <- matrix(data = 1, nrow = n, ncol = p + 1)
	for (i in 1:p) G[, p + 1 - i] <- ts[i:(i + n - 1)]
	lmo <- lmrob(formula = y ~ 0 + G, method = method, singular.ok = singular.ok, init = init)
	be <- lmo$coefficients
	for (i in 1:p) names(be)[i] <- paste("phi", i)
	names(be)[p + 1] <- "interc"
	return(list(coefficients = be, model = lmo))
}
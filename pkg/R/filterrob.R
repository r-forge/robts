filterrob <- function(x, filter = NULL, method = c("direct", "ar"), p = 0, na.action = na.fail, psi.l=2, psi.0=3) {
	method <- match.arg(method)
	if (!any(method == c("direct", "ar"))) stop("No valid method given.")
	x <- na.action(x)
	
	if (method == "direct") {
		if (!any(p == 1:length(x))) stop("No valid order of AR model given.")
		tss <- ARfilter(timeseries = x, p = p, psi.l=psi.l,psi.0=psi.0)[[5]][,p]
		resid <- NULL
	}
	
	if (method == "ar") {
		if (!is.numeric(filter)) stop("No valid vector of AR coefficients given.")
		res <- robfilterAR(timeseries = x, phi = filter, psi.l=psi.l, psi.0=psi.0)
		tss <- res$filtered.ts
		resid <- res$residuals
	}
	
	tss <- ts(tss)
	return(list(ts = tss, residuals = resid))
}

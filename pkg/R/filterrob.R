filterrob <- function(x, filter = NULL, method = c("direct", "ar"), p = 0, psifunc = smoothpsi, na.action = na.fail) {
	method <- match.arg(method)
	if (!any(method == c("direct", "ar"))) stop("No valid method given.")
	x <- na.action(x)
	
	if (method == "direct") {
		if (!any(p == 1:length(x))) stop("No valid order of AR model given.")
		tss <- ts.robfilter(ts = x, p = p, psifunc = psifunc)
		resid <- NULL
	}
	
	if (method == "ar") {
		if (!is.numeric(filter)) stop("No valid vector of AR coefficients given.")
		res <- robfilterAR(timeseries = x, phi = filter, psifunc = psifunc)
		tss <- res$filtered.ts
		resid <- res$residuals
	}
	
	tss <- ts(tss)
	return(list(ts = tss, residuals = resid))
}
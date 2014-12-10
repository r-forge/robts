## ARopt.filter - optimally fitted AR model by robust filtering
## Input:
##		ts: time series
##		pmax: maximum AR order considered
##		psifunc: a robust Psi-Function
## Output: order of optimal AR model

ARopt.filter <- function(tss, pmax = NULL, psifunc = smoothpsi) {
	tmax <- length(tss)
	if (!is.null(pmax)) if (pmax >= floor((tmax - 1) / 2)) {
		warning("Too less data for chosen pmax, corrected to greatest possible value.")
		pmax <- floor((tmax - 1) / 2) - 1
	}
	if (is.null(pmax)) pmax <- floor(min((tmax - 1) / 4, 10 * log(tmax, base = 10)))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	p <- pmax
	repeat {
		fits <- ARfilter(timeseries = tss, p = p, psifunc = psifunc)
		if (is.list(fits)) {
			RAICs <- fits[[5]]
			break
		} else p <- p - 1
		if (p == 0) stop("No successful computation for any p.")
	}
	popt <- which.min(RAICs)
	return(popt)
}
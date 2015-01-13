## robspectrum - Robust spectrum
## Input:
##		x: time series



robspectrum <- function(x, method = c("pgram", "ar"), psifunc = smoothpsi, spans = 8,
            acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"), 
            var1 = FALSE, sdestim = c("Qn", "scaleTau2", "mad", "sd"), plot =TRUE, ...) {
            	
	method <- match.arg(method)
	acf.fun <- match.arg(acf.fun)
	
	if (acf.fun == "acfrobfil") stop("This acf calculation function can not be used.")
	
	if (method == "pgram") {
		res <- robspec(x, psifunc = psifunc, acf.fun = acf.fun, spans = spans)
	}
	
	if (method == "ar") {
		acorf <- as.numeric(suppressWarnings(acfrob(x, fun = acf.fun, ...))$acf)
		res <- robspec.acf(acorf, tmax = length(x))
		if (var1) {sdestim <- function(y) {return(1)}} else {
			sdestim <- match.arg(sdestim)
			sdestim <- get(sdestim)
		}
		res$spec <- res$spec * sdestim(x)^2
	}
	
	res$method <- method
	class(res) <- "spec"
	if (plot) {
		plot(res)
		return(invisible(res))
	} else return(res)
}
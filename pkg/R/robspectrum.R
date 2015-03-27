## robspectrum - Robust spectrum
## Input:
##		x: time series



robspectrum <- function(x, method = c("pgram", "ar"), psifunc = smoothpsi, truncation = log(length(x),10)*10,bandwidth=FALSE,
            acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim", "nonrob"),
            arrob.method = c("yule-walker", "durbin-levinson", "robustregression", "filter", "gm", "nonrob"), 
            var1 = FALSE, sdestim = c("Qn", "scaleTau2", "mad", "sd"), plot =TRUE,kernel=c("parzen","daniell","bartlett","rectangular"), ...) {
	n <- length(x)
	kernel <- match.arg(kernel)
	if (kernel=="bartlett")
		fa <- 3/2
	if (kernel=="daniell")
    		fa <- 1
	if (kernel=="parzen")
		fa <- 1/0.539
	if (kernel=="rectangular")
		fa <- 1/2
        if (is.logical(bandwidth)) {
		bandwidth <- 2*pi/truncation*fa
		} else {if (bandwidth==0) truncation <- n-1 else {
				truncation <- 2*pi/bandwidth*fa
			}
		}	
		
	method <- match.arg(method)
	acf.fun <- match.arg(acf.fun)
	arrob.method <- match.arg(arrob.method)
	
	if (acf.fun == "acfrobfil") stop("This acf calculation function can not be used.")
	
	if (method == "pgram") {
		res <- robspec(x, psifunc = psifunc, acf.fun = acf.fun, truncation = truncation, arrob.method = arrob.method,
		kernel = kernel, smoothing = (bandwidth != 0))
	}
	
	if (method == "ar") {
		if (acf.fun == "nonrob") acorf <- as.numeric(suppressWarnings(acf(x, plot = FALSE))$acf[-1]) else
		acorf <- as.numeric(suppressWarnings(acfrob(x, fun = acf.fun, lag.max=truncation,..., plot = FALSE))$acf)
		res <- robspec.acf(acorf, tmax = length(x),M=truncation,kernel)
		if (var1) {sdestim <- function(y) {return(1)}} else {
			sdestim <- match.arg(sdestim)
			sdestim <- get(sdestim)
		}
		res$spec <- res$spec * sdestim(x)^2
	}
	
	res$method <- method
	res$bandwidth <- bandwidth
	class(res) <- "spec"
	if (plot) {
		plot(res)
		return(invisible(res))
	} else return(res)
}

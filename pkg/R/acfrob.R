## acfrob: robust acf calculation, wrapper function calling the acf...'s
## Arguments:
##		x: numeric vector: univariate time series
##		lag.max: numeric value of maximum lag, defaults to 10*log10(length(x))
##		fun: function used for computation
##		p: order of AR fit used if fun == "acfrobfil"
##		robfiltype: type of robust filtering if fun == "acfrobfil"
##		...: further arguments passed to the respective acf-function
## Return value: object of class acf

acfrob <- function(x, lag.max = NULL,
					fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA",
							"acfrank", "acfrobfil", "acftrim"),
					plot = TRUE, na.action = na.fail, p = NULL,
					robfiltype = c("emp", "pacf", "pacfMott"), ...) {
	
	fun <- match.arg(fun)
	robfiltype <- match.arg(robfiltype)
	
	# protective measures:
	
	x <- na.action(x)
	
	if (!any(fun == c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA",
						"acfrank", "acfrobfil", "acftrim"))) {
			warning("No valid acf function.")
			return(NA)
		}

	n <- length(x)
	
	if (is.null(lag.max)) lag.max <- 10 * log(n, base = 10)
	
	if (n < 2 + lag.max) {
		warning("You can't calculate this amount of lags. We will compute the maximal possible lag of n-2")
		lag.max <- n - 2
	}
	
	if (n < 4 * lag.max) {
		warning("Brockwell and Davis suggest to calculate only lags less n/4. Nevertheless we will calculate all lags you want, but you should be aware that the estimated acf for higher lags could be unreasonable.")
	}
	
	
	if (fun == "acfrobfil") {
		if (is.null(p)) {
			warning("No AR order given.")
			return(NA)
		}
		if (!any(robfiltype == c("emp", "pacf", "pacfMott"))) {
			warning("No valid type of robust filtering given.")
			return(NA)
		}
		acfval <- acfrobfil(timeseries = x, maxlag = lag.max, p = p, ...)[[which(robfiltype == c("emp", "pacf", "pacfMott"))]]
	} else {
		fun <- get(fun)
		acfval <- fun(timeseries = x, maxlag = lag.max, ...)
	}
	
	
	res <- list(lag = array(data = 1:lag.max, dim = c(lag.max, 1, 1)),
				acf = array(data = acfval, dim = c(lag.max, 1, 1)),
				type = "correlation",
				n.used = n,
				series = NULL,
				snames = NULL
				)
	class(res) <- "acf"
	if (plot) {
		plot(res)
		return(invisible(res))
	} else return(res)
}
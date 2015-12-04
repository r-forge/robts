## acfrob: robust acf calculation, wrapper function calling the acf...'s
## Arguments:
##		x: numeric vector: univariate time series
##		lag.max: numeric value of maximum lag, defaults to 10*log10(length(x))
##		fun: function used for computation
##		p: order of AR fit used if fun == "acfrobfil"
##		robfiltype: type of robust filtering if fun == "acfrobfil"
##		posdef: should autocorrelation be enforced to be positive semidefinite?
##		...: further arguments passed to the respective (p)acf-function
## Return value: object of class acf

acfrob <- function(x, lag.max = NULL,
					type = c("correlation", "partial"),
					fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA",
							"acfrank", "acfrobfil", "acftrim"),
					plot = TRUE, na.action = na.fail, p = NULL,
					robfiltype = c("emp", "pacf", "pacfMott"),
					partialtype = c("ranks", "durbin-levinson", "filter"),
					posdef=TRUE,...) {
	
	fun <- match.arg(fun)
	funname <- fun
	robfiltype <- match.arg(robfiltype)
	type <- match.arg(type)
	partialtype <- match.arg(partialtype)
	
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
	
	
	if (type == "correlation" | (type == "partial" & partialtype == "durbin-levinson")) {
	
		if (fun == "acfrobfil") {
			if (is.null(p)) {
				warning("No AR order given.")
				return(NA)
			}
			if (!any(robfiltype == c("emp", "pacf", "pacfMott"))) {
				warning("No valid type of robust filtering given.")
				return(NA)
			}
			acorf <- acfrobfil(timeseries = x, maxlag = lag.max, p = p, ...)
		} else {
			fun <- get(fun)
			acorf <- fun(timeseries = x, maxlag = lag.max, ...)
		}
		
		if (type == "partial") {
			acorf <- DurbinLev(acorf)$phis[[lag.max]]
		}
		
	}
	
	if (type == "partial" & partialtype != "durbin-levinson") {
		if (partialtype == "ranks") {
			acorf <- pacf2(x, lag.max = lag.max, ...)
		}
		if (partialtype == "filter") {
			acorf <- ARfilter(x, p = lag.max, ...)[[1]]
		}
	}
	
	if (type=="correlation") {
		if (posdef) {
			acorf2 <- try(acfposmaker(acorf),silent=TRUE)
			if (inherits(acorf2,"try-error")){
				warning("Transformation to a positiv definit acf failed.")
				} else {acorf <- acorf2}
		}
	}
	
	res <- list(lag = array(data = 1:lag.max, dim = c(lag.max, 1, 1)),
				acf = array(data = acorf, dim = c(lag.max, 1, 1)),
				type = type,
				n.used = n,
				series = deparse(substitute(x)),
				snames = NULL
				)
	class(res) <- "acf"
	if (plot) {
		confint <- konfband(fun=funname,n=n,...)
		plot(res,ci=0)
		if (length(confint==2)) {
			abline(h=confint[1],col="blue",lty="dashed")
			abline(h=confint[2],col="blue",lty="dashed")
			}
		return(invisible(res))
	} else return(res)
}

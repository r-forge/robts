## acfrob: robust acf calculation, wrapper function calling the acf...'s
## Arguments:
##		x: numeric vector: univariate time series
##		lag.max: numeric value of maximum lag, defaults to 10*log10(length(x))
##		approach: function used for computation
##		p: order of AR fit used if approach == "acfrobfil"
##		robfiltype: type of robust filtering if approach == "acfrobfil"
##		posdef: should autocorrelation be enforced to be positive semidefinite?
##		...: further arguments passed to the respective (p)acf-function
## Return value: object of class acf

acfrob <- function(x, lag.max = NULL,
  type = c("correlation", "partial"),
  approach = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA",
  "acfrank", "acfrobfil", "acftrim","acfregression"),
  plot = TRUE, na.action = na.fail,
  p = NULL, posdef = TRUE, ...) {
				#	robfiltype = c("filtered","ar"),
				#	partialtype = c("ranks", "durbin-levinson", "filter"),
  n <- length(x)
  
  # protective measures:
  type <- match.arg(type)	
	approach <- match.arg(approach)
  if (is.null(lag.max)) lag.max <- floor(10 * log(n, base = 10))
	x <- na.action(x) # handling of missing values
	if (n < 2 + lag.max) {
		warning("It is only possible to compute estimations lags up to n-2. Larger lags are not considered")
		lag.max <- n - 2
	}	
	if (n < 4 * lag.max) {
		warning("Brockwell and Davis (2006) suggest to calculate only lags less n/4. Nevertheless, we will calculate all lags you want but you should be aware that the estimated acf for higher lags could be unreasonable.")
	}
  
	funname <- approach
		
	if (type == "correlation" | (type == "partial" & partialtype == "durbin-levinson")) {
	
		if (approach == "acfrobfil") {
			if (is.null(p)) {
				p <- lag.max
			}
			if (!any(robfiltype == c("filtered", "ar"))) {
				warning("No valid type of robust filtering given.")
				return(NA)
			}
			acorf <- acfrobfil(timeseries = x, maxlag = lag.max, p = p, ...)
		} else {
			approach <- get(approach)
			acorf <- approach(timeseries = x, maxlag = lag.max, ...)
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
	
	# output:
	acfout <- list(lag = array(data = 1:lag.max, dim = c(lag.max, 1, 1)),
				acf = array(data = acorf, dim = c(lag.max, 1, 1)),
				type = type,
				n.used = n,
				series = deparse(substitute(x)),
				snames = NULL
				)
	class(acfout) <- "acf"
	if (plot) {
		confint <- konfband(approach=funname, n=n,...)
		plot(acfout, ci=0)
		if (length(confint==2)) {
			abline(h=confint[1], col="blue", lty="dashed")
			abline(h=confint[2], col="blue", lty="dashed")
		}
		return(invisible(acfout))
	} else {
    return(res)
  }
}

## acfrob: robust acf calculation, wrapper function calling the acf...'s
## Return value: object of class acf

acfrob <- function(x, lag.max = NULL,
  type = c("correlation", "covariance", "partial"),
  approach = c("GK", "median", "multi", "partrank", "RA",
  "rank", "filter", "trim", "bireg"), ...,
  plot = TRUE, na.action = na.fail, psd = TRUE, scalefn = Qn,
  partial.method = c("automatic", "durbin-levinson")) {
  
  # checks and preparations:
  type <- match.arg(type)	
	approach <- match.arg(approach)
	acffn <- get(paste("acfrob", approach, sep="."))
	partial.method <- match.arg(partial.method)
	series <- deparse(substitute(x))
	x <- na.action(as.ts(x)) # handling of missing values
	if (!is.numeric(x)) stop("'x' must be numeric")
  n <- length(x)
  if (is.null(lag.max)) lag.max <- floor(10 * log(n, base = 10))
  if (is.na(lag.max) || lag.max < 0) stop("'lag.max' must be at least 0")
	if (n < 2 + lag.max) {
		warning("It is only possible to compute estimations lags up to n-2. Larger lags are not considered")
		lag.max <- n - 2
	}	
	if (n < 4 * lag.max) {
		warning("Brockwell and Davis (2006) suggest to calculate only lags less n/4. Nevertheless, we will calculate all lags you want but you should be aware that the estimated acf for higher lags could be unreasonable.")
	}
  
  # actual calculations:
  if (type == "partial") {	
  	if (approach %in% c("filter", "partrank") & partial.method == "automatic") {
      pacfvalues <- acffn(x, lag.max = lag.max, partial = TRUE, ...)
    } else {
      acfvalues <- acffn(x, lag.max = lag.max, ...)
      pacfvalues <- DurbinLev(acfvalues)$phis[[lag.max]]
    }
    acfvalues <- pacfvalues    
  }	
  
	if (type %in% c("correlation", "covariance")) {
    acfvalues <- acffn(x, lag.max = lag.max, ...)
		if (psd) {
			acfvalues_psd <- try(make_acf_psd(acfvalues), silent=TRUE)
			if (inherits(acfvalues_psd, "try-error")){
				warning("Transformation to a positiv semidefinite acf failed. The returned result is not positive semidefinite.")
			} else {
        acfvalues <- acfvalues_psd
      }
		}
	}
	
	if (type == "covariance") {
    scaleest <- scalefn(x)
    acfvalues <- scaleest^2*acfvalues
  }
	
	# output:
	acf.out <- list(
    acf = array(data = c(1, acfvalues), dim = c(lag.max+1, 1, 1)),
    type = type,
    n.used = n,
    lag = array(data = 0:lag.max, dim = c(lag.max+1, 1, 1)),
    series = series,
    snames = NULL
	)
	class(acf.out) <- "acf"
	if (plot) {
		confint <- confband(approach=approach, n=n, ...)
		plot(acf.out, ci=0)
		if (length(confint==2)) {
			abline(h=confint[1], col="blue", lty="dashed")
			abline(h=confint[2], col="blue", lty="dashed")
		}
		return(invisible(acf.out))
	} else {
    return(acf.out)
  }
}

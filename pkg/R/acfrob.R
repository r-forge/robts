## acfrob: robust acf calculation, wrapper function calling the acf...'s
## Return value: object of class acf

acfrob <- function(x, lag.max = NULL,
  type = c("correlation", "covariance", "partial"),
  approach = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA",
  "acfrank", "acfrobfil", "acftrim", "acfregression"),
  plot = TRUE, na.action = na.fail, psd = TRUE, scalefn = Qn,
  partial.method = c("automatic", "durbin-levinson"), ...) {
  n <- length(x)
  
  # protective measures:
  type <- match.arg(type)	
	approach <- match.arg(approach)
	partial.method <- match.arg(partial.method)
  if (is.null(lag.max)) lag.max <- floor(10 * log(n, base = 10))
	x <- na.action(x) # handling of missing values
	if (n < 2 + lag.max) {
		warning("It is only possible to compute estimations lags up to n-2. Larger lags are not considered")
		lag.max <- n - 2
	}	
	if (n < 4 * lag.max) {
		warning("Brockwell and Davis (2006) suggest to calculate only lags less n/4. Nevertheless, we will calculate all lags you want but you should be aware that the estimated acf for higher lags could be unreasonable.")
	}
  
	acffn <- get(approach)

  if (type == "partial") {	
  	if (approach %in% c("acfrobfil", "acfpartrank") & partial.method == "automatic") {
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
    lag = array(data = 1:lag.max, dim = c(lag.max, 1, 1)),
    acf = array(data = acorf, dim = c(lag.max, 1, 1)),
    type = type,
    n.used = n,
    series = deparse(substitute(x)),
    snames = NULL
	)
	class(acf.out) <- "acf"
	if (plot) {
		confint <- confband(approach=funname, n=n,...)
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

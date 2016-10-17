arrob.filter <- function(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) {2*p}, na.action = na.fail, series = deparse(substitute(x)), psi.l = 2, psi.0 = 3) {
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
  x <- na.action(as.ts(x))
  if (anyNA(x)) stop("NAs in 'x'")
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (ists) xtsp <- tsp(x)
  xfreq <- frequency(x)
  x <- as.vector(x)
  n.used <- length(x)
  if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10)))) 
  if (order.max < 0L) stop("'order.max' must be >= 0")
  if (order.max >= n.used) stop("'order.max' must be < 'n.used'")
  cl <- match.call()
	n <- length(x)
	if (!is.null(order.max)) if (order.max >= floor((n - 1) / 2)) {
		warning("Not enough data for chosen model order 'order.max'. The largest possible value of 'order.max' is used.")
		order.max <- floor((n - 1) / 2) - 1
	}
	if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10))))
	if (order.max < 1) stop("Model order 'order.max' must be greater than zero. Try for example order.max = 1.")
	p_maxvalid <- order.max
	# Find the largest model order for which ARfilter can be applied:
	repeat {
		fits <- ARfilter(timeseries = x, p = p_maxvalid, aicpenalty=aicpenalty, psi.l=psi.l, psi.0=psi.0)
		if (is.list(fits)) {
			RAICs <- fits[[4]]
			break
		} else {
      p_maxvalid <- p_maxvalid - 1
    }
    if (p_maxvalid < order.max) warning("Could not fit models of order 'order.max'. Used largest value of p for which the model could be fitted instead.")
		if (p_maxvalid == 0) stop("No successful computation for any model order greater than zero.")
	}
	locandscale <- scaleTau2(x, mu.too=TRUE) # null model with p=0
	RAICs <- c(log(locandscale[2]^2)+aicpenalty(1)/n, RAICs, rep(NA, order.max-p_maxvalid)) # include the null model
	names(RAICs) <- 0L:order.max
	if(aic){
  	p_opt <- which.min(RAICs)[[1]] - 1	
	} else {
     p_opt <- p_maxvalid 
  }
 	if (p_opt==0) {
    coeff <- NULL
    var.pred <- locandscale[2]^2
    x.mean <- locandscale[1]
    resid <- x-x.mean
    partialacf <- rep(0, p_opt)
 	} else {
    coeff <- fits[[6]][p_opt,1:p_opt]
    var.pred <- fits[[2]][p_opt]^2
    x.mean <- scaleTau2(fits[[5]][, p_opt], mu.too=TRUE)[1]
    partialacf <- fits[[1]]
    xfilteredcen <- fits[[5]][, p_opt]-x.mean
    xfilteredcenma <- matrix(ncol=p_opt ,nrow=n-p_opt)
    for (i in 1:p_opt) {xfilteredcenma[, i] <- xfilteredcen[(p_opt-i+1):(n-i)]}
    resid <- c(rep(NA, p_opt), as.numeric(x[(p_opt+1):n]-x.mean-xfilteredcenma%*%coeff))	
 	}
  if (ists) {
        attr(resid, "tsp") <- xtsp
        attr(resid, "class") <- "ts"
    }
  res <- list(
		order = p_opt,
		ar = coeff,
		var.pred = var.pred,
		x.mean = x.mean,
		x.intercept = NULL,
		aic = RAICs, #the function ar returns the difference of the AIC values with the lowest one 
		n.used = n,
		order.max = p_maxvalid, #may be lower than the argument order.max in some cases 
		partialacf = array(partialacf, dim=c(length(partialacf), 1, 1)),
		resid = resid,
		method = "filter",
		series = series,
		frequency = xfreq,
		call = cl,
		asy.var.coef = NULL
	)
	class(res) <- "ar"
	return(res)
}

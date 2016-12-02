arrob.filter <- function(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p, na.action = na.fail, series = deparse(substitute(x)), psi.l = 2, psi.0 = 3) {
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
  x <- na.action(as.ts(x))
  if (anyNA(x)) stop("NAs in 'x'")
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (ists) xtsp <- tsp(x)
  xfreq <- frequency(x)
  x <- as.vector(x)
 	n <- length(x)
  if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10)))) 
  if (order.max < 0L) stop("'order.max' must be >= 0")
  if (order.max >= n) stop("'order.max' must be lower than the number of observations")
  cl <- match.call()
	if (!is.null(order.max)) if (order.max >= floor((n - 1) / 2)) {
		warning("Not enough data for chosen model order 'order.max'. The largest possible value of 'order.max' is used.")
		order.max <- floor((n - 1) / 2) - 1
	}
	if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10))))
	if (order.max < 1) stop("Model order 'order.max' must be greater than zero.")
	fits <- ARfilter(x, order.max=order.max, aicpenalty=aicpenalty, psi.l=psi.l, psi.0=psi.0)
	RAICs <- fits$aic # includes the null model
	names(RAICs) <- 0L:order.max
	order_selected <- if(aic) which.min(RAICs)[[1]] - 1 else max(which(!is.na(RAICs)))-1	
  coeff <- fits$ar[[order_selected+1]]
  var.pred <- fits$var[order_selected+1]
  x.mean <- scaleTau2(fits$filtered[, order_selected+1], mu.too=TRUE)[1]
  partialacf <- fits$pacf
  xfilteredcen <- matrix(ncol=order_selected+1, nrow=n-order_selected)
  for (i in 0:order_selected) {
    xfilteredcen[, i+1] <- (fits$filtered[, order_selected+1]-x.mean)[(order_selected-i+1):(n-i)]
  }
  resid <- c(rep(NA, order_selected), as.numeric(x[(order_selected+1):n] - x.mean - if(order_selected>0){xfilteredcen[, -1, drop=FALSE]%*%coeff}else{0}))
  if (ists) {
    attr(resid, "tsp") <- xtsp
    attr(resid, "class") <- "ts"
  }
  res <- list(
		order = order_selected,
		ar = coeff,
		var.pred = var.pred,
		x.mean = x.mean,
		x.intercept = NULL,
		aic = RAICs, #the function ar returns the difference of the AIC values with the lowest one 
		n.used = n,
		order.max = order.max,
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

arrob.yw <- function(x, order.max, aic = TRUE, aicpenalty = function(p) 2*p, na.action = na.fail, series = deparse(substitute(x)), acf.approach = c("GK", "median", "multi", "partrank", "RA", "rank", "filter", "trim", "bireg"), locfn = median, scalefn = Qn, ...) {
  cl <- match.call()
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (ists) xtsp <- tsp(x)
  xfreq <- frequency(x)
  x.original <- x
  x <- handle_missings_ts(x, na.action)
  n <- length(x) # actual number of observations used for estimation
  if (missing(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10)))) 
  if (order.max < 0L) stop("'order.max' must be >= 0")
  if (order.max >= n) stop("Argument 'order.max' must be lower than the length of the time series")
	if (order.max >= floor((n - 1) / 2)) {
		warning("Not enough data for chosen model order 'order.max'. The largest possible value of 'order.max' is used.")
		order.max <- floor((n - 1) / 2) - 1
	}
	acf.approach <- match.arg(acf.approach)
	
  RAICs <- rep(NA, order.max+1)
  names(RAICs) <- 0L:order.max
	xacf <- as.numeric(acfrob(x, lag.max = order.max, approach = acf.approach, plot = FALSE, ...)$acf)
	x.mean <- locfn(x)

	# null model:
	order_selected <- 0
  resid <- x - x.mean
	var.pred <- scalefn(resid)^2
	RAICs[1] <- log(var.pred) + aicpenalty(1)/n
  resid_selected <- resid  
	RAIC_selected <- RAICs[1] 
	coeff <- NULL
	partialacf <- rep(0, order.max)
  orders <- seq_len(order.max)
  ARparams <- acf2AR(xacf) 
	for (p in orders) {
		D <- matrix(nrow = n - p, ncol = p)
		for (j in 1L:p) D[, j] <- x[(p + 1 - j):(n - j)]
		D <- cbind(1, D)
		ph <- as.vector(ARparams[p, 1:p])
		ph1 <- c(x.mean * (1 - sum(ph)), ph)
		resid <- as.vector(x[(p + 1):n] - D %*% ph1)
		var.new <- scalefn(resid)^2
		RAICs[p+1] <- log(var.new) + aicpenalty(p+1)/(n-p)
		if (RAICs[p+1] < RAIC_selected || (!aic && p==order.max)) {
		  order_selected <- p
      RAIC_selected <- RAICs[p+1]
			resid_selected <- c(rep(NA, n-length(resid)), resid)
			var.pred <- var.new
			coeff <- ph
			partialacf <- ARMAacf(ar = coeff, lag.max = order.max, pacf = TRUE)
		}
	}
	resid_output <- naresid(attr(x, "na.action"), resid_selected)
  if (ists) {
        attr(resid_output, "tsp") <- xtsp
        attr(resid_output, "class") <- "ts"
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
		resid = resid_output,
		method = "Yule-Walker",
		series = series,
		frequency = xfreq,
		call = cl,
		asy.var.coef = NULL,
		x = x.original
	)
	attr(res, "na.action") <- attr(x, "na.action")
	class(res) <- c("arrob", "ar")
	return(res)
}

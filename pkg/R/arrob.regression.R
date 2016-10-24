arrob.regression <- function(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p, na.action = na.fail, series = deparse(substitute(x)), intercept = TRUE, scalefn = Qn, ...) {
  cl <- match.call()
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
  x <- na.action(as.ts(x))
  if (anyNA(x)) stop("NAs in 'x'")
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (ists) xtsp <- tsp(x)
  xfreq <- frequency(x)
  #x <- as.vector(x)
  n <- length(x)
  if (is.null(order.max)) order.max <- floor(min(c((n - 1) / 4, 10 * log(n, base = 10)))) 
  if (order.max < 0L) stop("'order.max' must be >= 0")
  if (order.max >= n) stop("Argument 'order.max' must be lower than the length of the time series")
	if (!is.null(order.max)) if (order.max >= floor((n - 1) / 2)) {
		warning("Not enough data for chosen model order 'order.max'. The largest possible value of 'order.max' is used.")
		order.max <- floor((n - 1 - as.numeric(intercept)) / 2) - 1
	}
	if (is.null(order.max)) order.max <- floor(min(c((n - 1 - as.numeric(intercept)) / 4, 10 * log(n, base = 10))))
	if (order.max < 1) stop("Model order 'order.max' must be greater than zero.")
	
  RAICs <- rep(NA, order.max+1)
  names(RAICs) <- 0L:order.max
  
	# null model:
	if (intercept) {
  	fit_selected <- lmrob(x ~ 1)#, ...)
  	x.intercept <- as.vector(fit_selected$coefficients)
  	x.mean <- x.intercept
    resid <- as.vector(fit_selected$residuals)
  	var.pred <- fit_selected$scale^2
  	RAICs[1] <- log(var.pred) + aicpenalty(1)/n
	} else {
    x.intercept <- NULL
    x.mean <- 0
    resid <- x
    var.pred <- scalefn(x)^2
    RAICs[1] <- log(var.pred)
  }
 	order_selected <- 0
 	RAIC_selected <- RAICs[1]
 	coeff <- NULL
 	partialacf <- rep(0, order.max)
  orders <- seq(along=numeric(order.max))
	for (p in orders) {
    y <- x[-(1:p)]
    G <- matrix(nrow = n-p, ncol = p)
    for (i in 1:p) G[, p + 1 - i] <- x[i:(i + n - p - 1)]
    if (intercept) {
      fit <- suppressWarnings(lmrob(y ~ G))#, ...))
    }	else {
      fit <- suppressWarnings(lmrob(y ~ 0 + G))#, ...))
    }
    if (fit$converged || p==order.max) {
      RAICs[p+1] <- log(fit$scale^2) + aicpenalty(p+as.numeric(intercept))/(n-p)
			if (RAICs[p+1] < RAIC_selected || (!aic && p==order.max)) {
				RAIC_selected <- RAICs[p+1]
				fit_selected <- fit
				order_selected <- p
				if (intercept) {
          coeff <- as.vector(fit$coefficients[1+(1:p)])
          x.intercept <- fit$coefficients[[1]]
          x.mean <- x.intercept/(1-sum(coeff))
        }	else {
          coeff <- as.vector(fit$coefficients[1:p])
        }
				var.pred <- fit$scale^2
				resid_selected <- c(rep(NA, p), as.vector(fit$residuals))
				partialacf <- ARMAacf(ar=coeff, lag.max=order.max, pacf=TRUE)
			}
    }
  }    
  if (ists) {
        attr(resid_selected, "tsp") <- xtsp
        attr(resid_selected, "class") <- "ts"
  }
  res <- list(
		order = order_selected,
		ar = coeff,
		var.pred = var.pred,
		x.mean = x.mean,
		x.intercept = x.intercept,
		aic = RAICs, #the function ar returns the difference of the AIC values with the lowest one 
		n.used = n,
		order.max = order.max,
		partialacf = array(partialacf, dim=c(length(partialacf), 1, 1)),
		resid = resid_selected,
		method = "regression",
		series = series,
		frequency = xfreq,
		call = cl,
		asy.var.coef = NULL
	)
	class(res) <- "ar"
	return(res)
}

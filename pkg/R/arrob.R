arrob <- function(x, aic = TRUE, order.max = NULL, method = c("yw", "regression", "gm", "filter"), na.action = na.fail, series = deparse(substitute(x)), ...) {
	method <- match.arg(method)
	res <- switch(method,
	 yw = arrob.yw(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 regression = arrob.regression(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 gm = arrob.gm(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 filter = arrob.filter(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...)
  )
  res$call <- match.call()	 
	return(res)
}

residuals.ar <- function(object, ...) object$resid
coef.ar <- function(object, ...) object$ar

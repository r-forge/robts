arrob <- function(x, aic = TRUE, order.max, method = c("yw", "regression", "gm", "filter"), na.action = na.fail, series = deparse(substitute(x)), ...) {

  # checks and preparations:
	method <- match.arg(method)
	# more checks are done by the function of the respective method

	res <- switch(method,
	 yw = arrob.yw(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 regression = arrob.regression(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 gm = arrob.gm(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 filter = arrob.filter(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...)
  )
  res$call <- match.call()	 
	return(res)
}

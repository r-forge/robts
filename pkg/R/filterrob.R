filterrob <- function(x, ar, ma, var.pred, method = c("statespace", "recursive", "statespaceARMA"), locfn, psi.l = 2, psi.0 = 3, na.action = na.fail, ...) {
  method <- match.arg(method)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
	
	if (missing(ar)) {
		arrobout <- arrob.filter(x, na.action = na.action, ...)
    ar <- arrobout$ar
		var.pred <- arrobout$var.pred
	} else {
    if(missing(var.pred)) stop("Argument 'var.pred' need to be provided. For a fitted AR model 'object' the residual variance be obtained by 'object$var.pred'.")
  }
	
	res <- switch(method,
    statespace = filterrob.statespace(x, ar = ar, var.pred = var.pred, locfn = locfn, psi.l = psi.l, psi.0 = psi.0, na.action = na.action),
    recursive = filterrob.recursive(x, ar = ar, var.pred = var.pred, locfn = locfn, psi.l = psi.l, psi.0 = psi.0, na.action = na.action),
    statespaceARMA = filterrob.statespaceARMA(x, ar = ar, ma = ma, var.pred = var.pred, locfn = locfn, psi.l = psi.l, psi.0 = psi.0, na.action = na.action)
  )
	return(res)
}

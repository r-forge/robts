filterrob <- function(x, ar = NULL, method = c("fit", "given"), psifn, locfn, scalefn, na.action = na.fail, ...) {
  method <- match.arg(method)
  if (!is.null(dim(x))) stop("Only implemented for univariate series")
  if(missing(psifn)) psifn <- function(x) M_psi(x, type="smooth")
	
	if (method == "fit") {
		ar <- arrob.filtered(x, na.action = na.action, ...)$ar
	}
	
	res <- filterrob.given(x, ar = ar, psifn = psifn, locfn = locfn, scalefn = scalefn, na.action = na.action)
	return(res)
}

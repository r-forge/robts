## ARopt.filter - optimally fitted AR model by robust filtering
## Input:
##		ts: time series
##		pmax: maximum AR order considered
##		psifunc: a robust Psi-Function
## Output: order of optimal AR model

ARopt.filter <- function(tss, pmax = NULL,aicpenalty=function(p) {return(2*p)},psi.l=2,psi.0=3) {
	tmax <- length(tss)
	if (!is.null(pmax)) if (pmax >= floor((tmax - 1) / 2)) {
		warning("Too less data for chosen pmax, corrected to greatest possible value.")
		pmax <- floor((tmax - 1) / 2) - 1
	}
	if (is.null(pmax)) pmax <- floor(min(c((tmax - 1) / 4, 10 * log(tmax, base = 10))))
	if (pmax < 1) stop("Too less data for reasonable model comparison. Try p = 1.")
	p <- pmax
	repeat {
		fits <- ARfilter(timeseries = tss, p = p, aicpenalty=aicpenalty,psi.l=psi.l,psi.0=psi.0)
		if (is.list(fits)) {
			RAICs <- fits[[4]]
			break
		} else p <- p - 1
		if (p == 0) stop("No successful computation for any p.")
	}
	# null model:
	locandscale <- scaleTau2(tss,mu.too=TRUE)
	RAICs <- c(log(locandscale[2]^2)+aicpenalty(1)/tmax, RAICs)
	popt <- which.min(RAICs)[1] - 1
	if (popt==0) {	coeff <- NULL
			var.pred <- locandscale[2]^2
			x.mean <- locandscale[1]
			resid <- tss-x.mean
			partialacf <- rep(0,pmax)
			} else {	coeff <- fits[[6]][popt,1:popt]
					var.pred <- fits[[2]][popt]^2
					x.mean <- scaleTau2(fits[[5]][,popt],mu.too=TRUE)[1]
					partialacf <- fits[[1]]
					tssfilteredcen <- fits[[5]][,popt]-x.mean
					tssfilteredcenma <- matrix(ncol=popt,nrow=tmax-popt)
					for (i in 1:popt) {tssfilteredcenma[,i] <- tssfilteredcen[(popt-i+1):(tmax-i)]}
					resid <- as.numeric(tss[(popt+1):tmax]-x.mean-tssfilteredcenma%*%coeff)	
				}
	return(list(order = popt, aic = RAICs,coefficients=coeff,var.pred=var.pred,x.mean=x.mean,resid=resid,partialacf=partialacf))
}

spectrumrob.pgram <- function(x, psifn = function(x) M_psi(x, type="smooth"), arrob.method = c("yw", "regression", "filter", "gm", "nonrobust"), kernel = c("parzen", "bartlett", "rectangular", "daniell"), truncation=round(log(length(x), 10) * 10), bandwidth = NULL, na.action = na.fail, series = deparse(substitute(x)), ...) {  
  
  # checks and preparations:
  arrob.method <- match.arg(arrob.method)
  kernel <- match.arg(kernel)
  series <- deparse(substitute(x))
  x <- handle_missings_ts(x, na.action)
	n <- length(x)
	
	#translate given bandwith to truncation or vice versa:
  fa <- switch(kernel, bartlett=3/2, daniell=1/0.90282336, parzen=1/0.539, rectangular=1/2)
  if(is.null(bandwidth)) {
    bandwidth <- 2*pi/truncation*fa  
  } else {
    truncation <- ifelse(bandwidth==0, n-1, round(2*pi/bandwidth*fa))
  }
  	
	if (arrob.method == "nonrobust") {
		arfit <- ar(x)
		resi <- as.numeric(matrix(arfit$resid[-(1:arfit$order)], ncol = 1))
		ka <- 1
	} else {
		arfit <- arrob(x, method = arrob.method, ...)
		resi <- psifn(as.numeric(arfit$resid[-(1:arfit$order)]) / sqrt(as.numeric(arfit$var.pred))) * sqrt(as.numeric(arfit$var.pred))
		inte <- function(x) psifn(x)^2 * dnorm(x)
		ka <- integrate(inte, -Inf, +Inf)$value
	}
	resi <- spec.taper(resi, p = 0.1)
	tmax <- length(resi)
	ph <- arfit$ar
	p <- length(ph)
	kmax <- floor(tmax / 2)
	ff <- frequency(x)
	ff <- seq.int(from = ff / n, by = ff / n, length.out = floor(n / 2))
	nff <- length(ff)
	XXre <- numeric(nff)
	XXim <- XXre
	sumRe <- XXre
	sumIm <- XXre
	for (k in seq(along=ff)) {
		XXre[k] <- sum(resi * cos(-2 * pi * ff[k] * 1:tmax))
		XXim[k] <- sum(resi * sin(-2 * pi * ff[k] * 1:tmax))
		sumRe[k] <- sum(ph * cos(2 * pi * ff[k] * 1:p))
		sumIm[k] <- sum(ph * sin(2 * pi * ff[k] * 1:p))
	}
	per <- (XXre^2 + XXim^2) / tmax
	if (bandwidth==0) { #no smoothing
    perS <- per
  } else {
		perS <- numeric(nff)
		spel <- c(rev(per),0,per,rev(per))
		weight <- do.call(kernel, args=list(k=c(-rev(ff),0,ff)*2*pi, window.type="spectrum", M=truncation))*2*pi/n
		for (i in 1:length(ff)) {
			perS[i] <- sum(weight*spel[i:(2*length(ff)+i)])
		}
	}
	Den <- (1 - sumRe)^2 + sumIm^2
	S <- perS / Den / ka
	res <- list(
    freq = ff,
    spec = S,
    coh = NULL,
    phase = NULL,
    series = series,
    snames = NULL,
    method = "pgramrob",
    bandwidth = bandwidth
  )
	attr(res, "na.action") <- attr(x, "na.action")
  class(res) <- "spec"
	return(res)
}

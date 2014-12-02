## robspec - Robust spectrum
## Input:
##		ts: time series
##		Arguments of smooth.fourier
## 		Arguments of lmrobARopt
## 		psifunc: Argument of ts.robfilter

robspec <- function(ts,
					smooth.type = "runmean", M = 1, taper = "cosine",
					method = "MM", singular.ok = FALSE, init = NULL, 
					psifunc = smoothpsi) {
	tsc <- ts - mean(ts)
	fitopt <- lmrobARopt(ts = tsc, interc = FALSE, method = method, singular.ok = singular.ok, init = init)
	p <- fitopt$p.optimal
	resfilt <- ts.robfilter(ts = fitopt$model$residuals, p = p, psifunc = psifunc)
	ph <- fitopt$coefficients
	per <- smooth.fourier(y = resfilt, smooth.type = smooth.type, M = M, taper = taper)
	tmax <- length(resfilt)
	kmax <- floor(tmax / 2)
	ff <- (-kmax):kmax / tmax
	sumRe <- numeric(2 * kmax + 1)
	sumIm <- sumRe
	for (k in seq_along(ff)) {
		sumRe[k] <- sum(ph * cos(2 * pi * ff[k] * 1:p))
		sumIm[k] <- sum(ph * sin(2 * pi * ff[k] * 1:p))
	}
	Den <- (1 - sumRe)^2 + sumIm^2
	S <- per / Den
	return(S)
}
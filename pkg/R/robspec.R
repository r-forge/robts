## robspec - Robust spectrum
## Input:
##		ts: time series
##		Arguments of smooth.fourier
## 		Arguments of lmrobARopt
## 		psifunc: Argument of ts.robfilter

robspec <- function(tss, psifunc = smoothpsi, acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"), spans = 8, 
				arrob.method = c("yule-walker", "durbin-levinson", "robustregression", "filter", "gm")) {
	acf.fun <- match.arg(acf.fun)
	arrob.method <- match.arg(arrob.method)
	stopifnot(is.numeric(tss), sum(is.na(tss)) == 0, spans %% 2 == 0)
	tmax <- length(tss)
	arfit <- arrob(x = tss, method = arrob.method, acf.fun = acf.fun)
	resi <- psifunc(arfit$resid / sqrt(arfit$var.pred))*sqrt(arfit$var.pred)
	resi <- spec.taper(resi, p = 0.1)
	tmax <- length(resi)
	ph <- arfit$ar
	p <- length(ph)
	kmax <- floor(tmax / 2)
	ff <- spectrum(tss, plot = FALSE)$freq
	XXre <- numeric(length(ff))
	XXim <- XXre
	sumRe <- XXre
	sumIm <- XXre
	for (k in seq_along(ff)) {
		XXre[k] <- sum(resi * cos(-2 * pi * ff[k] * 1:tmax))
		XXim[k] <- sum(resi * sin(-2 * pi * ff[k] * 1:tmax))
		sumRe[k] <- sum(ph * cos(2 * pi * ff[k] * 1:p))
		sumIm[k] <- sum(ph * sin(2 * pi * ff[k] * 1:p))
	}
	per <- (XXre^2 + XXim^2) / tmax
	if (spans == 0) perS <- per else {
		perS <- numeric(length(ff))
		for (i in seq_along(perS)) {
			mi <- max(i - spans / 2, 1)
			ma <- min(i + spans / 2, length(ff))
			ww <- c(1 / 2, rep(1, ma - mi - 1), 1 / 2)
			perS[i] <- weighted.mean(x = per[mi:ma], w = ww)
		}
	}
	Den <- (1 - sumRe)^2 + sumIm^2
	S <- perS / Den
	return(list(freq = ff, spec = S, coh = NULL, phase = NULL, series = NULL, snames = NULL, method = "AR"))
}

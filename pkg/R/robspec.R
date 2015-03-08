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
	ph <- arrob(x = tss, method = arrob.method, acf.fun = acf.fun)$ar
	p <- length(ph)
	D <- matrix(nrow = tmax - p, ncol = p)
	for (j in 1:p) D[, j] <- tss[(p + 1 - j):(tmax - j)]
	ph <- c(median(tss) * (1 - sum(ph)), ph)
	D <- cbind(1, D)
	resi <- tss[(p + 1):tmax] - D %*% ph
	resi <- psifunc(resi / Qn(resi))
	resi <- spec.taper(resi, p = 0.1)
	tmax <- length(resi)
	ph <- ph[-1]
	kmax <- floor(tmax / 2)
	ff <- 0:kmax / tmax
	XXre <- numeric(kmax + 1)
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
	perS <- numeric(kmax + 1)
	for (i in seq_along(perS)) {
		mi <- max(i - spans / 2, 1)
		ma <- min(i + spans / 2, kmax + 1)
		ww <- c(1 / 2, rep(1, ma - mi - 1), 1 / 2)
		perS[i] <- weighted.mean(x = per[mi:ma], w = ww)
	}
	Den <- (1 - sumRe)^2 + sumIm^2
	S <- perS / Den
	return(list(freq = ff, spec = S, coh = NULL, phase = NULL, series = NULL, snames = NULL, method = "AR"))
}
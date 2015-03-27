## robspec - Robust spectrum
## Input:
##		ts: time series
##		Arguments of smooth.fourier
## 		Arguments of lmrobARopt
## 		psifunc: Argument of ts.robfilter

robspec <- function(tss, psifunc = smoothpsi, acf.fun, truncation, 
				arrob.method, kernel, smoothing) {
	stopifnot(is.numeric(tss), sum(is.na(tss)) == 0)
	N <- length(tss)
	if (arrob.method == "nonrob" | acf.fun == "nonrob") {
		arfit <- ar(x = tss)
		arfit$resid <- matrix(arfit$resid[-(1:arfit$order)], ncol = 1)
	} else
	arfit <- arrob(x = tss, method = arrob.method, acf.fun = acf.fun)
	resi <- psifunc(as.numeric(arfit$resid) / sqrt(as.numeric(arfit$var.pred))) * sqrt(as.numeric(arfit$var.pred))
	resi <- spec.taper(resi, p = 0.1)
	tmax <- length(resi)
	ph <- arfit$ar
	p <- length(ph)
	kmax <- floor(tmax / 2)
	ff <- frequency(tss)
	ff <- seq.int(from = ff / N, by = ff / N, length.out = floor(N / 2))
	Nff <- length(ff)
	XXre <- numeric(Nff)
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
	if (!smoothing) perS <- per else {
		perS <- numeric(Nff)
		if (!any(kernel==c("parzen","bartlett","rectangular","daniell"))) {
			warning("This kernel is not implemented, using Parzen kernel instead.")
			kernel <- "parzen"
			}
		spel <- c(rev(per),0,per,rev(per))
		if (kernel=="parzen")
			gew <- parzen(c(-rev(ff),0,ff)*2*pi,window.type="spectrum",M=truncation)*2*pi/N
		if (kernel=="bartlett")
			gew <- bartlett(c(-rev(ff),0,ff)*2*pi,window.type="spectrum",M=truncation)*2*pi/N
		if (kernel=="rectangular")
			gew <- rectangular(c(-rev(ff),0,ff)*2*pi,window.type="spectrum",M=truncation)*2*pi/N
		if (kernel=="daniell")
			gew <- daniell(c(-rev(ff),0,ff)*2*pi,window.type="spectrum",M=truncation)*2*pi/N
		for (i in 1:length(ff)) {
			perS[i] <- sum(gew*spel[i:(2*length(ff)+i)])
			}
		}
	Den <- (1 - sumRe)^2 + sumIm^2
	S <- perS / Den
	return(list(freq = ff, spec = S, coh = NULL, phase = NULL, series = NULL, snames = NULL, method = "AR"))
}

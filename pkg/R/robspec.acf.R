## robspec.acf - Robust spectral density from ACF
## Input: acf
##		tmax: length of timeseries
## Output: robust spectral density

robspec.acf <- function(acf, tmax) {
	kmax <- floor(tmax / 2)
	hmax <- length(acf)
	ff <- 1:kmax / tmax
	s <- numeric(kmax)
	for (k in seq_along(ff)) {
		s[k] <- sum(acf * cos(2 * pi * ff[k] * 1:hmax))
	}
	s <- 2 + 4 * s
	return(list(freq = ff, spec = s, coh = NULL, phase = NULL, series = NULL, snames = NULL, method = "ACF"))
}
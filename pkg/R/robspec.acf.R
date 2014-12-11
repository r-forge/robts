## robspec.acf - Robust spectral density from ACF
## Input: acf
##		tmax: length of timeseries
## Output: robust spectral density

robspec.acf <- function(acf, tmax = length(acf) + 2) {
	kmax <- floor(tmax / 2)
	hmax <- length(acf)
	ff <- 0:kmax / tmax
	s <- numeric(kmax + 1)
	for (k in seq_along(ff)) {
		s[k] <- sum(acf * cos(2 * pi * ff[k] * 1:hmax))
	}
	s <- 2 + 4 * s
	names(s) <- ff
	return(s)
}
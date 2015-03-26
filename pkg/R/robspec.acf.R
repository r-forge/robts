## robspec.acf - Robust spectral density from ACF
## Input: acf
##		tmax: length of timeseries
## Output: robust spectral density

robspec.acf <- function(acf, tmax,M,kernel) {
	kmax <- floor(tmax / 2)
	hmax <- length(acf)
	ff <- 1:kmax / tmax
	s <- numeric(kmax)
	if (!any(kernel==c("parzen","bartlett","rectangular"))) {
		warning("This kernel is not implemented, using Parzen kernel instead.")
		kernel <- "parzen"
		}
	if (kernel=="parzen")
	weight <- parzen(1:hmax,M,window.type="acf")
	if (kernel=="bartlett")
	weight <- bartlett(1:hmax,M,window.type="acf")
	if (kernel=="rectangular") 
	weight <- 1
	if (M>hmax) {
		warning("Truncation point of acf is larger than calculated maximal lag of acf. This value is used instead")
		M <- hmax
		}
	for (k in seq_along(ff)) {
		s[k] <- sum(acf * weight * cos(2 * pi * ff[k] * 1:hmax))
	}
	s <- 1 + 2 * s
	return(list(freq = ff, spec = s, coh = NULL, phase = NULL, series = NULL, snames = NULL, method = "ACF"))
}


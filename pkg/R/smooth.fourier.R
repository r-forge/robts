## auxiliary function of robspec:
## robust smoothed spectral density estimate
## input:
##		y: data vector
##		smooth.type
##		M: half window width
##		taper: type of data taper
## output: spectral density estimate of length of y (+1)


smooth.fourier <- function(y, smooth.type = "runmean", M = 1, taper = "cosine") {
	tmax <- length(y)
	if (taper == "cosine") {
		sw <- 2 * pi / tmax
		tap <- seq(-pi + sw / 2, pi - sw / 2, sw)
		tap <- cos(tap) + 1
	}
	if (taper == "none") tap <- rep(1, tmax)
	yt <- y * tap
	kmax <- floor(tmax / 2)
	ff <- (-kmax):kmax / tmax
	XXre <- numeric(2 * kmax + 1)
	XXim <- XXre
	for (k in seq_along(ff)) {
		XXre[k] <- sum(yt * cos(-2 * pi * ff[k] * 1:tmax))
		XXim[k] <- sum(yt * sin(-2 * pi * ff[k] * 1:tmax))
	}
	SS <- (XXre^2 + XXim^2) / tmax
	if (smooth.type == "runmean") {
		SSs <- numeric(2 * kmax + 1)
		for (i in seq_along(SSs)) {
			m <- min(M, i - 1, 2 * kmax + 1 - i)
			SSs[i] <- mean(SS[(i - m):(i + m)])
		}
	}
	if (smooth.type == "none") SSs <- SS
	return(SSs)
}
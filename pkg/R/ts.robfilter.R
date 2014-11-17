##################
# ts.robfilter - robust filtering of time series
# input:
# 	ts: timeseries without NAs as vector
# 	p: order of used AR-fit
# 	psifunc: psi function
# output: robustly filtered time series
##################

ts.robfilter <- function(ts, p, psifunc = super) {
	res <- yohai(timeseries = ts, p = p, psifunc = psifunc)[[3]]
	return(res)
}
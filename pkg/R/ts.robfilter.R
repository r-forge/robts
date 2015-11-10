##################
# ts.robfilter - robust filtering of time series
# input:
# 	ts: timeseries without NAs as vector
# 	p: order of used AR-fit
# 	psifunc: psi function
# output: robustly filtered time series
##################

ts.robfilter <- function(ts, p) {
	res <- ARfilter(timeseries = ts, p = p)[[5]][,p]
	return(res)
}

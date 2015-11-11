##################
# ts.robfilter - robust filtering of time series
# input:
# 	ts: timeseries without NAs as vector
# 	p: order of used AR-fit
# 	psifunc: psi function
# output: robustly filtered time series
##################

ts.robfilter <- function(ts, p,psi.l=psi.l,psi.0=psi.0) {
	res <- ARfilter(timeseries = ts, p = p,psi.l=psi.l,psi.0=psi.0)[[5]][,p]
	return(res)
}

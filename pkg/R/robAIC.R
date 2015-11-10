##################
# robAIC - AIC of an AR(p) model
# input:
# 	ts: timeseries without NAs as vector
# 	p: order of used AR-fit
# 	psifunc: psi function
# output: AIC
##################

robAIC <- function(ts, p) {
	res <- ARfilter(timeseries = ts, p = p)[[5]][p]
	return(res)
}

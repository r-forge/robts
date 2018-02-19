## cusum for structural change in location
## input:
# x: time series
# var.method: how should the long run variance be calculated
# 	possibilities: 	window: uses a running window, see asymcusum for details
#			acf: estimates the acf of the timeseries, see asymacf for details
#			acfextra: estimates the acf by first two autocorrelations and extrapolation, see extracf for details
# overlapping: only used for option window in var.method: should windows overlap?
## output:
# res: complete trajectory of test statistic

changerob.cusum <- function(x, var.method = c("window", "acf", "acfextra"), overlapping = TRUE,shiftcorrect = TRUE, borderN = 10, ...){

  N <- length(x)
  b <- cumsum(x)
  b2 <- b[N]-b
  meandiff <- (b/(1:N)-b2/((N-1):0))[-N]
  t2 <- as.numeric((1:(N-1)))*as.numeric(((N-1):1))*meandiff/N^(1.5)
  if (shiftcorrect) {
  	tau <- which.max(t2[borderN:(N-borderN)])
  	jumpheight <- b[tau]/tau-(b[N]-b[tau])/(N-tau)
  	x[(tau+1):N] <- x[(tau+1):N]+jumpheight
  	}
  
  var.method <- match.arg(var.method)
  if (var.method=="window") {
  	asy <- asymvar.window(x=x,overlapping=overlapping,obs="untransformed",...)[[1]]
  	}
  if (var.method=="acf") {
  	asy <- asymvar.acf(x=x,obs="untransformed",...)[[1]]	
  	}
  if (var.method=="acfextra") {	
  	asy <- asymvar.acfextra(x=x,obs="untransformed",...)[[1]]
  	}
  res <- t2/sqrt(asy)
  return(res)
}

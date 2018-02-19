## wilcoxon-type cusum for structural change in location
## input:
# x: time series
# var.method: how should the long run variance be estimated
# 	possibilities: 	window: uses a running window, see asymcusum for details
#			acf: estimates the acf of the time series, see asymacf for details
#			acfextra: estimates the acf by first two autocorrelations and extrapolation, see extracf for details
# overlapping: only used for option window in var.method: should windows overlap?
## output:
# res: complete trajectory of test statistic

changerob.wilcox <- function(x, var.method = c("window", "acf", "acfextra"), overlapping = TRUE,shiftcorrect = TRUE, borderN = 10, ...){
  N <- length(x)
  
  t2 <- .Call("wilcoxsukz",x)[-N]
  
  if (shiftcorrect) {
  	tau <- which.max(t2[borderN:(N-borderN)])
  	jumpheight <- median(x[1:tau])-median(x[(tau+1):N])
  	x[(tau+1):N] <- x[(tau+1):N]+jumpheight
  	}
  
  var.method <- match.arg(var.method)
  if (var.method=="window") {
  	asy <- asymvar.window(x=x,overlapping=overlapping,obs="ranks",...)[[1]]
  	}
  if (var.method=="acf") {
  	asy <- asymvar.acf(x=x,obs="ranks",...)[[1]]	
  	}
  if (var.method=="acfextra") {	
  	asy <- asymvar.acfextra(x=x,obs="ranks",...)[[1]]
  	}
  
  res <- t2/sqrt(asy)
  return(res)
}


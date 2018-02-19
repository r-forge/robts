## Hodges Lehmann test for structural change in location
## input:
# x: time series
# var.method: how should the long run variance be estimated
# 	possibilities: 	window: uses a running window, see asymcusum for details
#			acf: estimates the acf of the timeseries, see asymacf for details
#			acfextra: estimates the acf by first two autocorrelations and extrapolation, see extracf for details
# overlapping: only used for option window in var.method: should windows overlap?	
## output:
# res: complete trajectory of teststatistic

changerob.HL <- function(x, var.method = c("window", "acf", "acfextra"), overlapping = TRUE, shiftcorrect = TRUE, borderN = 10, ...){
  N <- length(x)
  var.method <- match.arg(var.method)
  
  	threedots <- list(...)
  	index1 <- which(names(threedots) %in% c("dd","cc","K","type","order.max","aic","momentp"))
  	if (length(index1)==0) threedots1 <- list() else{
  		threedots1 <- threedots[index1]
  		}
  	index2 <- which(names(threedots) %in% c("type2","adjust","kernelused"))
  	if (length(index2)==0) threedots2 <- list() else{
  		threedots2 <- threedots[index2]
  		}
  
  diffi <- outer(x,x,"-")
  index <- upper.tri(diag(N))
  zeilennummer <- matrix(1:N,ncol=N,nrow=N)
  spaltennummer <- t(zeilennummer)
  
  diffi <- diffi[index]
  zeilennummer <- zeilennummer[index]
  spaltennummer <- spaltennummer[index]
  
  diffi <- rbind(diffi,zeilennummer,spaltennummer)
  index <- order(diffi[1,])
  diffi <- diffi[,index]
  t2 <- .Call("meddiffneu",diffi[1,],diffi[2,],diffi[3,],x[-N])
  
  t2 <- N^(-3/2)*t2*as.numeric((1:(N-1)))*as.numeric(((N-1):1))
  
  tau <- which.max(t2[borderN:(N-borderN)])
  
  if (shiftcorrect) {
  	jumpheight <- meddiff(x[1:tau],x[(tau+1):N])
  	x[(tau+1):N] <- x[(tau+1):N]+jumpheight
  	}
  
  u0 <- do.call(densdiff,c(list(x[1:tau]),list(x[(tau+1):N]),threedots2))
  
  if (var.method=="window") {
  	asy <- do.call(asymvar.window,c(list(x=x),list(overlapping=overlapping),list(obs="ranks"),threedots1))[[1]]
  	}
  if (var.method=="acf") {
  	asy <- do.call(asymvar.acf,c(list(x=x,obs="ranks"),threedots1))[[1]]	
  	}
  if (var.method=="acfextra") {	
  	asy <- do.call(asymvar.acfextra,c(list(x=x,obs="ranks"),threedots1))[[1]]
  	}
  
  res <- t2/sqrt(asy)*u0
return(res)
}


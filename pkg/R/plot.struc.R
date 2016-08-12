plot.struc <- function(x, ...){
  plot(x$trajectory,main=x$method,type="l",xlab="Time-index",ylab="Test trajectory",...,sub=paste("estimated change point: ",x$estimate,"; pvalue: ", round(x$p.value,3),sep=""))
  abline(v=x$estimate,lty="dotted",lwd=2,col="darkgreen")
  n <- length(x$critical.value)
  for (i in 1:n) {
  	abline(h=x$critical.value[i],lty="dashed",col="red")
	}
}

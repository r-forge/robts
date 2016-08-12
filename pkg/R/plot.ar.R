plot.ar <- function(x, ask=TRUE,...){
  n <- x$n.used
  p <- x$order
  devAskNewPage(ask=ask) 
  nn <- length(x$aic)
  plot(0:(nn-1),x$aic,main="AIC criterion",type="l",xlab="AR-order",ylab="AIC criterion",...,sub=paste("estimated AR-order: ",x$order,sep=""))
  abline(v=p,lty="dotted",lwd=2,col="darkgreen")
  acfrob(x$resid[(p+1):n], main="ACF of Pearson residuals")
  plot(x$resid,main="Residuals over time",xlab="Time index",ylab="Residuals")
  abline(h=qnorm(0.975,sd=sqrt(x$var.pred)),lty="dashed")
  abline(h=qnorm(0.025,sd=sqrt(x$var.pred)),lty="dashed")
  qqnorm(x$resid,ylab="Residual Quantiles")
  qqline(x$resid)
}

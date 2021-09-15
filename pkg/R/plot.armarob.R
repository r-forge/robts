plot.armarob <- function(x, ask = TRUE, ci = 0.95, ...){
  n <- x$n.used
  p <- max(x$order[1],x$order[2]+1)
  devAskNewPage(ask = ask)   
  
  #Information criterion for choosing the model order:
  nn1 <- length(x$aic[,1])
  nn2 <- length(x$aic[1,])
  if ((nn1>0)&(nn2>0)) {
  image(x=0:(nn2-1),y=0:(nn1-1),z=x$aic,xlab="AR order", ylab="MA order",main = "Values of robust AIC criterion")
  }
  
  #Robust ACF of the residuals:
  acf_resid <- acfrob(x$resid[(p+1):n], plot = FALSE)
  plot(acf_resid, ci = ci, main = "Robust ACF of the residuals")
  
  #Residuals against time:
  plot(x$resid, main = "Residuals over time", xlab = "Time", ylab = "Residuals")
  if(ci > 0) abline(h = c(-1, +1) * qnorm((ci+1)/2, sd = sqrt(x$var.pred)), lty = "dashed", col="blue")
  
  #QQ plot:
  qqnorm(as.numeric(x$resid), ylab = "Normal QQ plot of the residuals")
  qqline(x$resid)
}

plot.ar <- function(x, ask = TRUE, ci = 0.95, ...){
  n <- x$n.used
  p <- x$order
  devAskNewPage(ask = ask)   
  
  #Information criterion for choosing the model order:
  nn <- length(x$aic)
  plot(0:(nn-1), x$aic, main = "Information criterion for different model orders", type = "l", xlab = "Model order", ylab = "Information criterion", sub = paste("(chosen order: ", x$order, ")", sep=""))
  abline(v = p, lty = "dashed", lwd = 2, col = "blue")
  
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

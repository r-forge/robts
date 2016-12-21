plot.changerob <- function(x, ...){
  plot(x$trajectory, main = x$method, type = "l", xlab = "Time", ylab = "Test statistic", ylim = range(x$trajectory, x$critical.value, na.rm = TRUE), ..., sub = paste("(estimated change point: ", x$estimate, "; pvalue: ", round(x$p.value,3), ")", sep=""))
  abline(v = x$estimate, lty = "solid", col = "red")
 	abline(h = x$critical.value, lty = "dashed", col = "blue")
}

acfregression <- function(x, lag.max, regression.method=c("MM", "lts", "quan")) {
  regression.method <- match.arg(regression.method)
  acfvalues <- numeric(lag.max)
  n <- length(x)
  for (i in 1:lag.max) {
  if (regression.method=="MM") acfvalues[i] <- coef(lmrob(x[(i+1):n]~x[1:(n-i)]))[2]
  if (regression.method=="lts") acfvalues[i] <- coef(ltsReg(x[(i+1):n]~x[1:(n-i)]))[2]
  if (regression.method=="quan") acfvalues[i] <- coef(rq(x[(i+1):n]~x[1:(n-i)]))[2]
  }
  return(acfvalues)
}




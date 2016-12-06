acfrob.bireg <- function(x, lag.max, regression.method = c("MM", "LTS", "L1")) {
  regression.method <- match.arg(regression.method)
  acfvalues <- numeric(lag.max)
  n <- length(x)
  for (i in 1:lag.max) {
    if (regression.method=="MM") acfvalues[i] <- coef(lmrob(x[(i+1):n] ~ x[1:(n-i)]))[2]
    if (regression.method=="LTS") acfvalues[i] <- coef(ltsReg(x[(i+1):n] ~ x[1:(n-i)]))[2]
    if (regression.method=="L1") acfvalues[i] <- coef(rq(x[(i+1):n] ~ x[1:(n-i)]))[2]
  }
 	res <- list(
   acfvalues = acfvalues,
   are = NA
  )
  	
  return(res)
}

acfregression <- function(timeseries,maxlag,regression.method="MM") {
x <- timeseries
erg <- numeric(maxlag)
n <- length(x)
for (i in 1:maxlag) {
if (regression.method=="MM") erg[i] <- coef(lmrob(x[(i+1):n]~x[1:(n-i)]))[2]
if (regression.method=="lts") erg[i] <- coef(ltsReg(x[(i+1):n]~x[1:(n-i)]))[2]
if (regression.method=="quan") erg[i] <- coef(rq(x[(i+1):n]~x[1:(n-i)]))[2]
}
return(erg)
}




###############
# auxililiary function: calculates a reweighted version of the Qn estimator
# input
# x: first variable to calculate correlation
# boundq: 
# output:
###############

reweightedQn <- function(x, boundq = 0.999) {
  n <- length(x)
  s0 <- Qn(x)
  m0 <- median(x)
  robmalahanobis <- ((x-m0)/s0)^2
  ordermalahanobis <- sort(robmalahanobis)
  
  bound <- qchisq(boundq, 1)	# bound from which distribution function of malahanobis distances is compared
  dn <- max((pchisq(ordermalahanobis, df = 1)-rank(ordermalahanobis)/n+1/n)*(ordermalahanobis>=bound)) # maximum positiv! difference between theoretical and empirical distribution function, compared from borderquantile (+1/n is result of jumpdiscontinuity of empirical distribution) 
  
  alphan <- n - floor(dn*n)		# calculation of the cutoff value
  alphan <- ordermalahanobis[alphan]
  weight <- (robmalahanobis < alphan)	 # cut all observations which are greater  then the observation where the difference between the distribution functions was maximal
  gooddata <- x[weight==1]
  
  result <- sd(gooddata)
  return(result)
}

BurgM <- function(ut, vt, eps=10^(-3), maxit=100) {
  n <- length(ut)
  sn <- (Qn(ut)^2+Qn(vt)^2)/2
  if (sn==0) {
  	warning("Variance estimation is 0")
  	return(NA)
  	}
  w <- function(x) 3/(1+x)
  fn <- corGK(ut,vt)
  fnal <- fn+1
  i <- 1
  while(abs(fn-fnal)>eps){
  fnal <- fn
  dn <- ut^2+vt^2-2*fn*ut*vt
  wn <- w(dn/sn)
  fn <- 2*sum(wn*ut*vt)/sum(wn*(ut^2+vt^2))
  if(fn==Inf){
  	warning("estimated variance is 0")
  	return(NA)
  	}
  sn <- 1/2/(n)*sum(wn*(ut^2+vt^2))
  i <- i+1
  if (i > maxit) {
  	warning("Iteration does not converge")
  	return(fn)
  	}
  }
  return(fn)
}
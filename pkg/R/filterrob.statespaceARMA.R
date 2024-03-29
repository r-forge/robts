filterrob.statespaceARMA <- function(x, ar, ma, var.pred, locfn, psi.l = 2, psi.0 = 3, na.action = na.fail) {
  ma <- -ma
  ists <- is.ts(x)
  if (ists) xtsp <- tsp(x)
  x <- handle_missings_ts(x, na.action)
  if (!is.numeric(x)) stop("'x' must be numeric")  
  if(missing(locfn)) locfn <- function(x) scaleTau2(x, mu.too = TRUE)[1]
  n <- length(x)
  p <- length(ar)
  q <- length(ma)
  kd <- max(p,q+1)
  mu <- locfn(x)
  x_centered <- x - mu
  if((p==0)&(kd==0)) {
    filtered <- M_psi(x_centered/var.pred, type = "smooth", k = c(psi.l, psi.0))
    residuals <- x_centered
  } else {
    # calculating first autocorrelation:
    k <- psi.l
    l <- psi.0
    a <- (2*k^2*l^2)/(k-l)^3
    b <- -(l^3+k*l^2+4*k^2*l)/(k-l)^3
    d <- (2*l^2+2*k*l+2*k^2)/(k-l)^3
    e <- -(l+k)/(k-l)^3
    
    # consistency correction for residual variance (like scaleTau2):
    conscorr <- function(c1=4.5, c2=3) {
      b <- c2*qnorm(3/4)
      corfa <- 2*((1-b^2)*pnorm(b)-b*dnorm(b)+b^2)-1
      return(corfa)
    }
   
    acvf <- tacvfARMA(phi=ar, theta=ma, maxLag = kd-1, sigma2 = var.pred)
    filterout <- .Call("filter4", c(x_centered, 0, rep(0,n)), ar, ma, var.pred, a, b, d, e, k, l, 4.5, 3, conscorr(), acvf)
    filtered <- filterout[1:n]
    residuals <- filterout[n+1+(1:n)]
    for(i in (q+1):n) {residuals[i] <- residuals[i]+t(ma)%*%residuals[(i-1):(i-q)]}
    residuals[1:q] <- NA
  }
  filtered <- filtered + mu
  residuals <- naresid(attr(x, "na.action"), residuals)
  filtered <- napredict(attr(x, "na.action"), filtered)
  if (ists) {
    attr(residuals, "tsp") <- attr(filtered, "tsp") <- xtsp
    attr(residuals, "class") <- attr(filtered, "class") <- "ts"
  } 
  res <- list(filtered = filtered, residuals = residuals)
 	attr(res, "na.action") <- attr(x, "na.action")
  return(res)
}

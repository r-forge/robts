## ARsolveYW - calculation of AR model coefficients by solving the Yule-Walker equations
## input:
## 	acf: autocorrelation function starting from lag 0
## 	p: order of the AR model
## output: AR model coefficients

ARsolveYW <- function(autocor, p) {
	p <- as.integer(p)
	stopifnot(p >= 1)
	#if (length(acf) < p) stop("The acf need to be available up to lag 'p'")
	r <- autocor[1+(1:p)]
	R <- diag(nrow = p)
	if (p > 1) {
    for (i in 1:(p-1)) {
  		R[i, (i+1):p] <- r[1:(p-i)]
  		R[p-i+1, 1:(p-i)] <- r[(p-i):1]
  	}
	}
	ph <- solve(a = R, b = r)  #
	return(ph)
}

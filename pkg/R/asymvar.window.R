## estimates long run variance (see Dehling et al. 2013) (overlapping and nonoverlapping block sums)
## input:
# x: time series
# dd: how should block length be calculated
#	possibilities:	"independent": block-length-factor independent of correlation structure of the sample
#			"carlstein-cor": block-length-factor dependent of acf(1) of transformed time series
#			"carlstein-Qn": blocklength-factor dependent of acf(1) of time series (robustly estimated)
# rhotrue: acf(1) of time series
# obs: if obs=untransformed, the long-run-variance for a cusum test is calculated, if obs="ranks" the long-run-variance for a wilcoxon changepoint test is calculated
# p: which centered mean should be used, see Peligrad and Shao (1995) for details
## output:
# l: used block length of variance estimator
# er: estimated long run variance

asymvar.window <- function(x, overlapping = TRUE, obs = c("untransformed", "ranks"), dd = c("independent", "carlstein-cor", "carlstein-Qn"), momentp = 1){
  N <- length(x)
  dd <- match.arg(dd)
  obs <- match.arg(obs)
  if (obs=="ranks") x=rank(x)/length(x)

  l <- round((3*N)**(1/3)+1)
 
  if (dd=="carlstein-cor"){rho=cor(x[1:(N-1)],x[2:N])
    l <- max(ceiling((N*(2*rho/(1-rho*rho))**2)**(1/3)),1)
  }
  if (dd=="carlstein-Qn"){rho=acfrob.GK(x, lag.max=1)
    l <- max(ceiling((N*(2*rho/(1-rho*rho))**2)**(1/3)),1)
  }
  if (sum(dd==c("carlstein-cor","carlstein-Qn","independent"))==0) {
    warning(paste(dd,"is not implemented. Using blocklength independent of correlation structure instead."))
  }
  phibar <- mean(x)
  if (overlapping) {
    S <- .Call("runmean",x,l)[1:(N-l)]
    S <- (abs(S-l*phibar)/sqrt(l))^momentp
    cp <- 2^(-momentp/2)*sqrt(pi)/gamma((momentp+1)/2)
    er <- sum(S)/(N-l)*cp
    res <- list(lrv=er^(2/momentp),blocklength=l)
  } else {
    k <- floor(N/l)
    xma <- matrix(x[1:(k*l)],ncol=k,nrow=l)
    S <- apply(xma,2,sum)
    S <- (abs(S-l*phibar)/sqrt(l))^momentp
    cp <- 2^(-momentp/2)*sqrt(pi)/gamma((momentp+1)/2)
    er <- sum(S)/k*cp
    res <- list(lrv=er^(2/momentp),blocklength=l)
  }
 return(res)
}

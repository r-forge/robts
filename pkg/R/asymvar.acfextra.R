## estimates long run variance by extrapolating the acf
## input:
# y: timeseries
# tau: time of structural break if timeseries is not stationary
# h: jump height
# order.max: number of estimated autocorreelations which are used for extrapolation (if aic=TRUE maximal possible AR order)
# aic: should AR order be estimated based on aic criterium?
## output:
# asy: estimated long run variance

asymvar.acfextra <-function(x,shiftcorrect=TRUE,obs=c("untransformed","ranks"),borderN=10,order.max=2,aic=FALSE,...){
 n=length(x)
 obs <- match.arg(obs)
 if (obs=="untransformed") {
 	if (shiftcorrect) {
		zz=rep(0,n)
		for (jj in borderN:(n-borderN)){
			zz[jj]=jj*(n-jj)*abs(mean(x[1:jj])-mean(x[(jj+1):n]))
			}
		tauh=which.max(zz==max(zz))
		height= mean(x[1:tauh])-mean(x[(tauh+1):n])
		x[(tauh+1):n] <- x[(tauh+1):n]+height
		}
	}
if (obs=="ranks") {
	if(shiftcorrect) {
		zz=rep(0,n)
		for (jj in borderN:(n-borderN)){
			zz[jj]=jj*(n-jj)*abs(meddiff(x[1:jj],x[(jj+1):n]))
			}
		tauh=which.max(zz)
		height= meddiff(x[1:tauh],x[(tauh+1):n])
		x[(tauh+1):n] = x[(tauh+1):n]+height
		}
  	x=edf(x)
	}
 armodel <- ar(x,order.max=order.max,aic=aic)
 if (length(armodel$ar)==0) acfest <- rep(0,n-1) else{
 	acfest <- ARMAacf(ar=armodel$ar,lag.max=n-1)[-1]
	}
 w=(n-1:(n-1))/n
 asy=1+2*sum(acfest*w)
 asy=asy*var(x)
 erg <- list(lrv=asy,order=length(armodel$ar))	
 return(erg)
}

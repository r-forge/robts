## estimates long run variance (see Dehling et al. 2013) (overlapping and nonoverlappingblocksums Si)
## input:
# x: timeseries
# dd: how should block length be calculated
#	possibilities:	"independent": blocklength-factor independent of correlation structure of the sample
#			"carlstein-cor": blocklength-factor dependent of acf(1) of transformed timeseries
#			"carlstein-Qn": blocklength-factor dependent of acf(1) of timeseries (robustly estimated)
# rhotrue: acf(1) of timeseries
# shiftcorrect: should variace estimation consider a possible levelshift
# obs: if obs=untransformed, the longrunvariance for a cusum test is calculated, if obs="ranks" the longrunvariance for a wilcoxon changepoint test is calculated
# borderN: if shiftcorrect=TRUE, how many observations at the beginning and end are left out for a change point search
# p: which centered mean should be used, see (Peligrad and Shao 1995 for details)
## output:
# l: used blocklength of variance estimator
# er: estimated long run variance

asymvar.window=function(x=x,overlapping=FALSE,shiftcorrect=TRUE,obs=c("untransformed","ranks"),dd=c("independent","carlstein-cor","carlstein-Qn"),borderN=10,momentp=1,...){
 N=length(x)
 dd <- match.arg(dd)
 obs <- match.arg(obs)
 if (obs=="untransformed") {
	if (shiftcorrect) {
		zz=rep(0,N)
		for (jj in borderN:(N-borderN)){
			zz[jj]=jj*(N-jj)*abs(mean(x[1:jj])-mean(x[(jj+1):N]))
			}
		tauh=which.max(zz==max(zz))
		height= mean(x[1:tauh])-mean(x[(tauh+1):N])
		x[(tauh+1):N] <- x[(tauh+1):N]+height
		}
	}
 if (obs=="ranks") {
	if(shiftcorrect) {
		zz=rep(0,N)
		for (jj in borderN:(N-borderN)){
			zz[jj]=jj*(N-jj)*abs(meddiff(x[1:jj],x[(jj+1):N]))
			}
		tauh=which.max(zz)
		height= meddiff(x[1:tauh],x[(tauh+1):N])
		x[(tauh+1):N] = x[(tauh+1):N]+height
		}
  	x=edf(x)
	}

 l=round((3*N)**(1/3)+1)
 
 if (dd=="carlstein-cor"){rho=cor(x[1:(N-1)],x[2:N])
   l=max(ceiling((N*(2*rho/(1-rho*rho))**2)**(1/3)),1)
 }
 if (dd=="carlstein-Qn"){rho=acfGK(x,maxlag=1)
   l=max(ceiling((N*(2*rho/(1-rho*rho))**2)**(1/3)),1)
 }
 if (sum(dd==c("carlstein-cor","carlstein-Qn","independent"))==0)
   {warning(paste(dd,"is not implemented. Using blocklength independent of correlation structure instead."))}

 phibar=mean(x)
 if (overlapping) {
 S=x
 for (i in 1:(N-l)){S[i]=sum(x[i:(i+l-1)])}
   S=(abs(S-l*phibar)/sqrt(l))^momentp
   cp=2^(-momentp/2)*sqrt(pi)/gamma((momentp+1)/2)
   er=sum(S)/(N-l)*cp
   erg <- list(lrv=er^(1/momentp),blocklength=l)
 return(erg)
 }
 if (!overlapping) {
 k=floor(N/l)
 S=rep(0,k)
 for (i in 1:k){S[i]=sum(x[((i-1)*l+1):(i*l)])}
 S=(abs(S-l*phibar)/sqrt(l))^momentp
 cp=2^(-momentp/2)*sqrt(pi)/gamma((momentp+1)/2)
 er=sum(S)/k*cp
 erg <- list(lrv=er^(1/momentp),blocklength=l)
 return(erg)
 }
 }

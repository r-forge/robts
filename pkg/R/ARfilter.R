#######################
# auxiliary function: calculates the partial correlations using robustly filtered values
# input:
# timeseries: timeseries without NAs as vector
# p: order of AR fit
# psifunc: psifunc: a robust Psi-Function
# output
# partial: partialautocorrelations
# varold: estimated variance of innovations
# xhat: robustly filtered timeseries
# phiacf: calculates acf until lag p (like proposed by oja (partial correlation with tue use of ranks)) 
########################

ARfilter <- function(timeseries,p,psifunc=smoothpsi,aicpenalty=function(n,p) {return(2*p/n)}) {
n <- length(timeseries)
aicv <- numeric(p)
partial <- numeric(p)		# partial autocorrelations
phiacf <- numeric(p)

timeseriesalt <- matrix(nrow=n,ncol=p)
varold <- numeric(p)
aralt <- matrix(ncol=p,nrow=p)
acfalt <- matrix(ncol=p,nrow=p)

# calculating first autocorrelation

segam <- function(x)  filterinit(x,timeseries,psifunc=psifunc)[[1]]
op <- try(optimize(segam,c(-1,1)),silent=TRUE)
if (inherits(op,"try-error")){
	warning("Optimization failed. Partial Autocorrelation can not be computed.")
	return(NA)
}
partial[1] <- op$minimum
helppar <- partial[1]
aralt[1,1] <- partial[1]
varold[1] <- op$objective
phiacf[1] <- partial[1]		# estimated acf
aicv[1] <- log(varold[1])+aicpenalty(n-1,1)
# calculation of the other autocorrelations
timeseriesalt[,1] <- filterinit(partial[1],timeseries,psifunc=psifunc)[[2]]


if (p > 1) for (j in 2:p) {
	segam <- function(x) filter(x,helppar,varold[j-1],timeseries,psifunc=psifunc)[[1]]
	op <- try(optimize(segam,c(-1,1)),silent=TRUE)
	if (inherits(op,"try-error")){
		warning("Optimization failed. Partial Autocorrelation can not be computed.")
		return(partial)
	}	
	partial[j] <- op$minimum
	timeseriesalt[,j] <- try(filter(partial[j],helppar,varold[j-1],timeseries,psifunc=psifunc)[[2]],silent=TRUE)
	if (inherits(timeseriesalt[,j],"try-error")){
		warning("Calculation failed, partial correlation can not be computed.")
		return(NA)
	}	
	varold[j] <- op$objective
	Phi <- numeric(j)
	Phi[j] <- partial[j]
	
	# updating the memory j predictors
	
	for (i in 1:(j-1)) {
		Phi[i] <- helppar[i]-partial[j]*helppar[j-i]
		}
	helppar <- Phi
	phiacf[j] <- sum(helppar[1:(j-1)]*phiacf[(j-1):1])+Phi[j]*(1-sum(helppar[1:(j-1)]*phiacf[1:(j-1)]))	# acf recursion like in the paper of oja (we use it also in function acfpartrank)
	aicv[j] <- log(varold[j])+aicpenalty(n-j,j)	

	aralt[j,1:j] <- helppar 
}

erg <- list(partial,varold,phiacf,aicv,timeseriesalt,aralt)
names(erg) <- c("partial autocorrelations","variance of innovations","autocorrelation","aic","robustly filtered timeseries","AR-coefficients")
return(erg)
}

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

yohai <- function(timeseries,p,psifunc=super) {
n <- length(timeseries)
aicv <- numeric(p)
partial <- numeric(p)		# partial autocorrelations
phiacf <- numeric(p)

timeseriesalt <- matrix(nrow=n,ncol=p)
aralt <- matrix(ncol=p,nrow=p)
acfalt <- matrix(ncol=p,nrow=p)

# calculating first autocorrelation

segam <- function(x)  sega(x,timeseries,psifunc=psifunc)[[1]]
op <- try(optimize(segam,c(-1,1)),silent=TRUE)
if (inherits(op,"try-error")){
	warning("Optimization failed. Partial Autocorrelation can not be computed.")
	return(NA)
}
partial[1] <- op$minimum
helppar <- partial[1]
varold <- op$objective
phiacf[1] <- partial[1]		# estimated acf
aicv[1] <- log(varold)+2/(n-1)
# calculation of the other autocorrelations
timeseriesalt[,1] <- sega(partial[1],timeseries,psifunc=psifunc)[[2]]


for (j in 2:p) {
	segam <- function(x) segaII(x,helppar,varold,timeseries,psifunc=psifunc)[[1]]
	op <- try(optimize(segam,c(-1,1)),silent=TRUE)
	if (inherits(op,"try-error")){
		warning("Optimization failed. Partial Autocorrelation can not be computed.")
		return(partial)
	}	
	partial[j] <- op$minimum
	varold <- op$objective
	Phi <- numeric(j)
	Phi[j] <- partial[j]
	
	# updating the memory j predictors
	
	for (i in 1:(j-1)) {
		Phi[i] <- helppar[i]-partial[j]*helppar[j-i]
		}
	helppar <- Phi
	phiacf[j] <- sum(helppar[1:(j-1)]*phiacf[(j-1):1])+Phi[j]*(1-sum(helppar[1:(j-1)]*phiacf[1:(j-1)]))	# acf recursion like in the paper of oja (we use it also in function Autokor)
aicv[j] <- log(varold)+2*j/(n-j)	
timeseriesalt[,j] <- try(segaII(0,partial[1:j],varold,timeseries,psifunc=psifunc)[[2]],silent=TRUE)
if (inherits(timeseriesalt[,j],"try-error")){
	warning("Calculation failed, partial correlation can not be computed.")
	return(NA)
}
aralt[1:j,j] <- helppar 
}

erg <- list(partial,varold,segaII(0,partial,varold,timeseries,psifunc=psifunc)[[2]],phiacf,aicv,timeseriesalt)
names(erg) <- c("partial autocorrelations","variance of innovations","robustly filtered timeseries","autocorrelation","aic","timeseries")
return(erg)
}
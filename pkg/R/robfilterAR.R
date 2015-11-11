#####################
# robfilterAR: calculates a robustly filtered timeseries
# input:
#	 timeseries: timeseries without NA as vector
#	 phi: coefficients of AR model fitted
#	 psifunc: a robust Psi-Function
# output: named list of 2:
#	filtered.ts: robustly filtered timeseries
#	residuals: residuals of AR fit
#####################


robfilterAR <- function(timeseries, phi,psi.l=psi.l,psi.0=psi.0) {
psifunc <- function(x) return(smoothpsi(x,k=psi.l,l=psi.0))

n <- length(timeseries)
p <- length(phi)

# state transition matrix
Phi <- matrix(data = 0, nrow = p, ncol = p)
Phi[1, ] <- phi
if (p > 1) diag(Phi[2:p, 1:(p - 1)]) <- 1

mu <- median(timeseries)
xhat <- numeric(n+p)		# robustly predicted values		
ug <- numeric(n)		# residuals
sigmau <- Qn(timeseries)^2

# Estimation of a recursion-start value for the filtering error covariance

datamatrix <- matrix(nrow=n-p+1,ncol=p)
for (i in 1:p) {
	datamatrix[,i] <- timeseries[(p-i+1):(n-i+1)] 
	}
P <- try(covMcd(datamatrix)$cov,silent=TRUE)
if (inherits(P,"try-error")) {
	warning("start estimation of covariance filtering errors failed.")
	return(NA)	
	}

# calculation of robustly filteres values

for(i in (p+1):(n+p)) {
	xhat[i] <- mu+rev(Phi[1,])%*%(xhat[(i-p):(i-1)]-mu)	# one step ahead predictions
	ug[i-p] <- timeseries[i-p]-xhat[i]			# residuals
M <- Phi%*%P%*%t(Phi)						# recursion for covariance of prediction error
M[1,1] <- M[1,1]+sigmau
mhat <- M[,1]
if (mhat[1]<=0) {
	warning("Prediction error variance is less or equal 0.")
	return(NA)
	}
shat <- mhat[1]^0.5
P <- M-1/shat^2*psifunc(ug[i-p]/shat)/(ug[i-p]/shat)*mhat%*%t(mhat)	# recursion for variance of filtering error
xhat[i:(i-p+1)] <- xhat[i:(i-p+1)]+mhat/shat*psifunc(ug[i-p]/shat)		# robustly filtered values
}

erg <- list(filtered.ts = xhat[-(1:p)], residuals = ug)
return(erg)
}

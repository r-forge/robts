#####################
# auxiliary function: calculates the robustly filtered timeseries and its robust variance (k-th step of Durbin-Levionson algorithm)
# input
# parauto: k-th partial autocorrelation
# oldmemory: memory k-1 predictors
# oldvar: variance of m-1 prediction
# timeseries: timeseries without NA as vector
# psifunc: a robust Psi-Function
# output: robustly filtered timeseries and its robust variance (k-th step)
#####################


segaII <- function(parauto,oldmemory,oldvar,timeseries,psifunc=super) {


# protective measures

if (abs(parauto>1)) {
	warning("Partial Autocorrelation can not have an absolute value larger then 1.")
}

p <- length(oldmemory)+1
n <- length(timeseries)

# updating the memory k predictors

Phi <- numeric(p)
Phi[p] <- parauto
for (i in 1:(p-1)) {
	Phi[i] <- oldmemory[i]-parauto*oldmemory[p-i]
	}
Phi <- rbind(Phi,cbind(diag(rep(1,p-1)),0))

mu <- median(timeseries)
xhat <- numeric(n+p)		# robustly predicted values		
ug <- numeric(n)		# residuals


sigmau <- oldvar^2*(1-Phi[1,p]^2)	# variance of memory m prediczion errors

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

sdhat <- scaleTau2(ug) 
erg <- list(sdhat,xhat[-(1:p)])
return(erg)
}

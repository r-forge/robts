#####################
# auxiliary function: calculates the robustly filtered timeseries and its robust variance (first step of Durbin-Levionson algorithm)
# input
# autocor1: first lag autocorrelation
# timeseries: timeseries without NA as vector
# psifunc: a robust Psi-Function
# output: robustly filtered timeseries and its robust variance (first step)
#####################

sega <- function(autocor1,timeseries,psifunc=super) {

# protective measure

if (abs(autocor1>1)) {
	warning("Partial Autocorrelation can not have an absolute value larger then 1.")
}

n <- length(timeseries)
mu <- median(timeseries)
xhat <- numeric(n+1)		# robustly predicted values
ug <- numeric(n)		# residuals
xhat[1] <- 0			# start of recursion

# variance estimation
sigma <- scaleTau2(timeseries)^2	
sigmau <- sigma*(1-autocor1^2)
P <- sigma			# start of recursion

# calculation of robustly filteres values

for(i in 2:(n+1)) {
xhat[i] <- mu+autocor1*(xhat[i-1]-mu)	# one step ahead predictions
ug[i-1] <- timeseries[i-1]-xhat[i]		# residuals
m <- autocor1^2*P+sigmau			# recursion for variance of prediction error
if (m<=0) {
	warning("Prediction error variance is less or equal 0.")
	return(NA)
	}
shat <- m^0.5
P <- m-psifunc(ug[i-1]/shat)/(ug[i-1]/shat)*m		# recursion for variance of filtering error
xhat[i] <- xhat[i]+shat*psifunc(ug[i-1]/shat)	# robustly filtered value
}

# variance estimation

sdrob <- scaleTau2(ug) 				
erg <- list(sdrob,xhat[-1])
return(erg)
}
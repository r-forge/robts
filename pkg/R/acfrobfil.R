##################
# estimation of the acf using robustly filtered values
# input
# timeseries: timeseries without NAs as vector
# p: order of used AR-fit
# maxlag: maximal lag of interest
# output
# acft: acf using the empirical acf of the robustly filtered timeseries
# acfp: acf using the estimated partial autocorrelations
# acfp2: acf using the estimated partial autocorrelations (transformation like in ojas paper robust autocorrelation using robust autocorrelation)
##################


acfrobfil <- function(timeseries,p,maxlag,psifunc=smoothpsi) {


n <- length(timeseries)



# estimating the partial autocorrelation and robustly filtered values

estimate <- try(ARfilter(timeseries,p,psifunc),silent=TRUE)
if (inherits(estimate,"try-error")){
	warning("Calculation of the acf failed.")
	return(NA)
	}
if (length(estimate)==1) {
warning("Calculation of the acf failed.")
	return(NA)
}

# estimating acf using robustly filtered values

robfiltered <- estimate[[6]][,p]
acft <- try(acf(robfiltered,plot=FALSE,lag.max=maxlag)$acf[-1],silent=TRUE)
if (inherits(acft,"try-error")){
	warning("Calculation of the acf failed.")
	return(NA)
	}



# estimate acf using partial autocorrelations

acfp2 <- numeric(maxlag+1)		# Oja-Type
acfp2[1] <- 1
acfp2[2:(p+1)] <- estimate[[4]][1:p]

acfp <- numeric(maxlag+1)		# Yule-Walke-Equiation-Type
acfp[1] <- 1
werte <- estimate[[1]][1:p]		# partial autocorrelations

if (p>1){
acfp[2:(p+1)] <- solveYuleII(werte)} else{acfp[2] <- werte}	# autocorrelations until lag p

# further autocorrelations if needed using yule walker recursion

if (maxlag > p) {
	for (i in (p+1):(maxlag)) {
		rho <- 0
		rho2 <- 0
		for (j in 1:p){
			rho <- rho+werte[j]*acfp[abs(i-j)+1]
			rho2 <- rho2+werte[j]*acfp2[abs(i-j)+1]
			} 
		acfp[i+1] <- rho
		acfp2[i+1] <- rho2
		}
	}

result <- list(acft,acfp[-1],acfp2[-1],p)

p <- which.min(estimate[[5]])
# estimating acf using robustly filtered values

robfiltered <- estimate[[6]][,p]
acftaic <- try(acf(robfiltered,plot=FALSE,lag.max=maxlag)$acf[-1],silent=TRUE)
if (inherits(acftaic,"try-error")){
	warning("Calculation of the acf failed.")
	return(NA)
	}


# estimate acf using partial autocorrelations

acfp2aic <- numeric(maxlag+1)		# Oja-Type
acfp2aic[1] <- 1
acfp2aic[2:(p+1)] <- estimate[[4]][1:p]

acfpaic <- numeric(maxlag+1)		# Yule-Walke-Equiation-Type
acfpaic[1] <- 1
werte <- estimate[[1]][1:p]		# partial autocorrelations

if (p>1){
acfpaic[2:(p+1)] <- solveYuleII(werte)} else{acfpaic[2] <- werte}	# autocorrelations until lag p

# further autocorrelations if needed using yule walker recursion

if (maxlag > p) {
	for (i in (p+1):(maxlag)) {
		rho <- 0
		rho2 <- 0
		for (j in 1:p){
			rho <- rho+werte[j]*acfpaic[abs(i-j)+1]
			rho2 <- rho2+werte[j]*acfp2aic[abs(i-j)+1]
			} 
		acfpaic[i+1] <- rho
		acfp2aic[i+1] <- rho2
		}
	}

result <- list(acft,acfp[-1],acfp2[-1],acftaic,acfpaic[-1],acfp2aic[-1],p,robfiltered)
return(result)
}
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


acfrobfil <- function(timeseries,p,maxlag,psifunc=smoothpsi,robfiltype="filtered",aic=TRUE) {


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

if (aic) p <- which.min(estimate[[4]])

# estimating acf using robustly filtered values

if (robfiltype=="filtered"){

	robfiltered <- estimate[[5]][,p]
	acfv <- try(acf(robfiltered,plot=FALSE,lag.max=maxlag)$acf[-1],silent=TRUE)
	if (inherits(acfv,"try-error")){
		warning("Calculation of the acf failed.")
		return(NA)
		}
	}

# estimate acf using AR fit

if (robfiltype=="ar") {
	acfv <- ARMAacf(ar=estimate[[6]][,p],lag.max=maxlag)[-1]
	}
return(acfv)
}

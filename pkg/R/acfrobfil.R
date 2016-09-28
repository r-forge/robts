##################
# estimation of the acf using robustly filtered values
# input
# x: time series without NAs as vector
# p: order of used AR-fit
# lag.max: maximal lag of interest
# output
# acft: acf using the empirical acf of the robustly filtered x
# acfp: acf using the estimated partial autocorrelations
# acfp2: acf using the estimated partial autocorrelations (transformation like in ojas paper robust autocorrelation using robust autocorrelation)
##################


acfrobfil <- function(x, lag.max, p = lag.max, robfil.method = c("filtered", "ar"), aic = TRUE, psi.l = 2, psi.0 = 3) {
  n <- length(x)
  robfil.method <- match.arg(robfil.method)
  
  # estimating the partial autocorrelation and robustly filtered values:
    estimate <- try(ARfilter(x, p, psi.l=psi.l, psi.0=psi.0), silent=TRUE)
  if (inherits(estimate, "try-error") || length(estimate)==1){
  	stop("Calculation of the acf failed.")
 	}
  
  if (aic) p <- which.min(estimate[[4]])
  
  #estimating acf using robustly filtered values:
  if (robfil.method=="filtered"){
 	  robfiltered <- estimate[[5]][, p]
  	acfvalues <- try(acf(robfiltered, plot=FALSE, lag.max=lag.max)$acf[-1], silent=TRUE)
  	if (inherits(acfvalues, "try-error")){
  		stop("Calculation of the acf failed.")
 		}
 	}
  
  # estimate acf using AR fit:
  if (robfil.method=="ar") {
  	acfvalues <- ARMAacf(ar=estimate[[6]][, p], lag.max=lag.max)[-1]
  	}
  	
  return(acfvalues)
}

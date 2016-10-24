##################
# estimation of the acf using robustly filtered values
# input
# x: time series without NAs as vector
# order.max: order of used AR-fit
# lag.max: maximal lag of interest
# output
# acft: acf using the empirical acf of the robustly filtered x
# acfp: acf using the estimated partial autocorrelations
# acfp2: acf using the estimated partial autocorrelations (transformation like in ojas paper robust autocorrelation using robust autocorrelation)
##################


acfrob.filter <- function(x, lag.max, order.max = lag.max, robfil.method = c("filtered", "ar"), aic = TRUE, aicpenalty=function(p) {2*p}, psi.l = 2, psi.0 = 3, partial = FALSE) {
  n <- length(x)
  robfil.method <- match.arg(robfil.method)
  
  # estimating the partial autocorrelation and robustly filtered values:
    estimate <- try(ARfilter(x, order.max=order.max, aicpenalty=aicpenalty, psi.l=psi.l, psi.0=psi.0), silent=TRUE)
  if (inherits(estimate, "try-error") || length(estimate)==1){
  	stop("Calculation of the acf failed.")
 	}
  
  order_selected <- if (aic) which.min(estimate$aic)-1 else order.max
  
  #estimating acf using robustly filtered values:
  if (robfil.method=="filtered"){
 	  robfiltered <- estimate$filtered[, order_selected+1]
  	acfvalues <- try(acf(robfiltered, plot=FALSE, lag.max=lag.max, type=ifelse(partial, "partial", "correlation"))$acf[-1], silent=TRUE)
  	if (inherits(acfvalues, "try-error")){
  		stop("Calculation of the (p)acf failed.")
 		}
 	}
  
  # estimate acf using AR fit:
  if (robfil.method=="ar") {
    if (partial) {
      pacfvalues <- estimate$pacf
      return(pacfvalues)
    } 
  	acfvalues <- if(order_selected > 0) ARMAacf(ar=estimate$ar[[order_selected+1]], lag.max=lag.max)[-1] else  rep(0, lag.max)
 	}
  	
  return(acfvalues)
}

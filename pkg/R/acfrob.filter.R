acfrob.filter <- function(x, lag.max, order.max = lag.max, robfil.method = c("filtered", "ar"), aic = TRUE, aicpenalty=function(p) {2*p}, psi.l = 2, psi.0 = 3, partial = FALSE) {
  n <- length(x)
  robfil.method <- match.arg(robfil.method)
  
  # estimating the partial autocorrelation and robustly filtered values:
    estimate <- suppressWarnings(try(ARfilter(x, order.max=order.max, aicpenalty=aicpenalty, psi.l=psi.l, psi.0=psi.0), silent=TRUE))
  if (inherits(estimate, "try-error") || length(estimate)==1){
  	stop("Filtering of the time series failed.")
 	}
  
  #selecting the order of the AR model:
  if (aic) {
    order_selected <- which.min(estimate$aic)-1
  } else {
    order_selected <- max(which(!is.na(estimate$aic)))-1
    if (order_selected < order.max) warning(paste("It was not possible to fit an AR model of the desired order ", order.max, ". Instead an\nAR model of the highest possible order ", order_selected, " was used for filtering the time series.", sep=""))
  }
  
  #estimating acf of the residuals using robustly filtered values:
  if (robfil.method=="filtered"){
 	  robfiltered <- estimate$filtered[, order_selected+1]
  	acfvalues <- try(acf(robfiltered, plot=FALSE, lag.max=lag.max, type=ifelse(partial, "partial", "correlation"))$acf[as.numeric(!partial)+(1:lag.max), 1, 1], silent=TRUE)
  	if (inherits(acfvalues, "try-error")){
  		stop("Calculation of the (p)acf of the filtered values failed.")
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
 	
 	res <- list(
   acfvalues = acfvalues,
   are = NA
  )
  	
  return(res)
}

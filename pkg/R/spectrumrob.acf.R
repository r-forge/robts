spectrumrob.acf <- function(x, acf.approach = c("GK", "median", "multi", "partrank", "RA", "rank", "trim", "bireg", "nonrobust"), kernel = c("parzen", "bartlett", "rectangular", "daniell"), truncation = log(length(x), 10) * 10, bandwidth = NULL, var1 = FALSE, scalefn = Qn, na.action = na.fail, series = deparse(substitute(x)), ...) {
  
  # checks and preparations:
  acf.approach <- match.arg(acf.approach)
  kernel <- match.arg(kernel)
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x)) # handling of missing values
  n <- length(x)

	#translate given bandwith to truncation or vice versa:
  fa <- switch(kernel, bartlett=3/2, daniell=1/0.90282336, parzen=1/0.539, rectangular=1/2)
  if(is.null(bandwidth)) {
    bandwidth <- 2*pi/truncation*fa  
  } else {
    truncation <- ifelse(bandwidth==0, n-1, round(2*pi/bandwidth*fa))
  }
	
	if (acf.approach == "nonrobust") {
    acorf <- as.numeric(suppressWarnings(acf(x, lag.max=truncation, plot = FALSE))$acf[-1])
  } else {
    acorf <- as.numeric(suppressWarnings(acfrob(x, lag.max=truncation, approach = acf.approach, ..., plot = FALSE))$acf[-1])
  }
   	
	hmax <- length(acorf)
	if(hmax != truncation){
    warning(paste("The autocorrelation function could not be calculated for all", truncation, "lags as given\nby the argument 'truncation' (or as implied by the argument 'bandwidth').\nInstead the maximal number of lags", hmax, "was used. The list element 'bandwidth'\nin the output reflects the actual bandwith corresponding with this value."))
    truncation <- hmax
    bandwidth <- 2*pi/truncation*fa
  }
  kmax <- floor(n / 2)  
 	ff <- 1:kmax / n
	s <- numeric(kmax)
	weight <- do.call(kernel, args=list(k=1:hmax, window.type="acf", M=truncation))
  for (k in seq(along=ff)) {
		s[k] <- sum(acorf * weight * cos(2 * pi * ff[k] * 1:hmax))
	}
	s <- 1 + 2 * s
	if (!var1) s <- s * scalefn(x)^2
	res <- list(
    freq = ff,
    spec = s,
    coh = NULL,
    phase = NULL,
    series = series,
    snames = NULL,
    method = "acfrob",
    bandwidth = bandwidth
  )
	class(res) <- "spec"  	
	return(res)		
}

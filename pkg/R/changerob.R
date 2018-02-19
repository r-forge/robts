### robust changepoint tests ###
## robust tests for a change in location or scale
# input:
# x: time series
# property: whether to test against a change in location or scale
# test: which test to use: Cusum test, Wilcoxon-test or Hodges-Lehmann test
# conf.level: desired significance level of the test
# alternative: power against which alternative
# 	possibilities:	two.sided: classical approach
#			increase: mean increases
#			decrease: mean decreases
# var.method: procedure to estimate the long run variance
#	possibilities:	window: subsampling procedure
#			acf: kernel estimator
#			acfextra: extrapolates acf from ar-fit
# overlapping: should distinct or overlapping blocks be used to calculate the long run variance, only used if var.method=window
# shiftcorrect: If TRUE the estimation of the lon run variance considers a location change
# plot: draw trajectory of test statistic

changerob <- function(x, property = c("location", "scale"), test = c("HL", "Wilcoxon", "CUSUM"), conf.level = 0.95, alternative = c("two.sided", "increase", "decrease"), var.method = c("window", "acf", "acfextra"), overlapping = TRUE, shiftcorrect = TRUE, borderN = 10, plot = FALSE, ...) {

  test <- match.arg(test)
  property <- match.arg(property)
  alternative <- match.arg(alternative)
  var.method <- match.arg(var.method)
  n <- length(x)
  if (shiftcorrect) {
  	if (borderN > n/4) warning(paste("Argument 'borderN' was too large and is set to ", floor(n/4), ".", sep=""))
  	borderN <- min(floor(n/4), borderN)
  }
  
  if (property=="scale") {
  	x <- abs(x-median(x))
  	x0 <- x==0
  	xmin <- min(x[!x0])
  	x[x0] <- xmin 
  	}
  
  if (test=="CUSUM") {
  	testtrajectory <- changerob.cusum(x=x, var.method=var.method, overlapping=overlapping,shiftcorrect=shiftcorrect, borderN=borderN, ...)
  }
  if (test=="Wilcoxon") {
  	testtrajectory <- changerob.wilcox(x=x, var.method=var.method, overlapping=overlapping,shiftcorrect=shiftcorrect, borderN=borderN, ...)
  }
  if (test=="HL") {
  	testtrajectory <- changerob.HL(x=x, var.method=var.method,overlapping=overlapping,shiftcorrect=shiftcorrect, borderN=borderN, ...)
  }
  
  data.name <- names(x)
  
  if (alternative=="two.sided") {
  	p.value <- 1-pKS(max(abs(testtrajectory)))
  	statistic <- max(abs(testtrajectory))
  	estimate <- which.max(abs(testtrajectory))
  	cval <- c(-1,+1)*qKS(conf.level)
 	}
  
  if (alternative=="decreae") {
  	p.value <- exp(-max(c(0,testtrajectory))^2)
  	statistic <- max(testtrajectory)
  	estimate <- which.max(testtrajectory)
  	cval <- sqrt(-log(1-conf.level))
 	}
  
  if (alternative=="increase") {
  	p.value <- exp(-min(c(0,testtrajectory))^2)
  	statistic <- min(testtrajectory)
  	estimate <- which.min(testtrajectory)
  	cval <- -sqrt(-log(1-conf.level))
 	}
  names(estimate) <- "estimated change point"
  names(statistic) <- "test statistic"
  null.value <- 0
  names(null.value) <- "possible level shift"
  res <- list(
    statistic = statistic,
    p.value = p.value,
    estimate = estimate,
    null.value = null.value,
    alternative = alternative,
    method = paste(test, "test for change in", property),
    data.name = data.name,
    trajectory = testtrajectory,
    critical.value = cval
  )
  attr(res$critical.value, "conf.level") <- conf.level
  class(res) <- c("changerob", "htest")
  if (plot) plot(res)
  return(res)
}


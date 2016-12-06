## wilcoxon-type cusum for structural change in location
## input:
# z: timeseries
# varmethod: how should the long run variance be calculated
# 	possibilities: 	window: uses a running window, see asymcusum for details
#			acf: estimates the acf of the timeseries, see asymacf for details
#			acfextra: estimates the acf by first two autocorrelations and extrapolation, see extracf for details
# overlapping: only used for option window in varmethod: should windows overlap?
## output:
# t2: complete trajectory of teststatistic

strucchange.wilcox <- function(x,varmethod=c("window","acf","acfextra"),overlapping=TRUE,shiftcorrect=TRUE,borderN=10,...){
N <- length(x)

t2 <- .Call("wilcoxsukz",x)[-N]

if (shiftcorrect) {
	tau <- which.max(t2[borderN:(N-borderN)])
	jumpheight <- median(x[1:tau])-median(x[(tau+1):N])
	x[(tau+1):N] <- x[(tau+1):N]+jumpheight
	}

varmethod <- match.arg(varmethod)
if (varmethod=="window") {
	asy <- asymvar.window(x=x,overlapping=overlapping,obs="ranks",...)[[1]]
	}
if (varmethod=="acf") {
	asy <- asymvar.acf(x=x,obs="ranks",...)[[1]]	
	}
if (varmethod=="acfextra") {	
	asy <- asymvar.acfextra(x=x,obs="ranks",...)[[1]]
	}

t2 <- t2/sqrt(asy)
return(t2)
}


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
varmethod <- match.arg(varmethod)
if (varmethod=="window") {
	asy <- asymvar.window(x=x,overlapping=overlapping,shiftcorrect=shiftcorrect,obs="ranks",borderN=borderN,...)[[1]]
	}
if (varmethod=="acf") {
	asy <- asymvar.acf(x=x,shiftcorrect=shiftcorrect,obs="ranks",...)[[1]]	
	}
if (varmethod=="acfextra") {	
	asy <- asymvar.acfextra(x=x,shiftcorrect=shiftcorrect,obs="ranks",...)[[1]]
	}
t2 <- numeric(N-1)
for (m in 1:(N-1)){
	x2=x[(m+1):N]
	x1=x[1:m]
	n=N-m
	sx <- wilcox.test(x1,x2)$statistic
	t0=2*sx-m*n
	t2[m]=sqrt(m*n)*t0/(sqrt(m+n)*N*N)	
	}
t2 <- t2/asy
return(t2)
}


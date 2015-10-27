## Hodges Lehmann test for structural change in location
## input:
# z: timeseries
# varmethod: how should the long run variance be calculated
# 	possibilities: 	window: uses a running window, see asymcusum for details
#			acf: estimates the acf of the timeseries, see asymacf for details
#			acfextra: estimates the acf by first two autocorrelations and extrapolation, see extracf for details
# overlapping: only used for option window in varmethod: should windows overlap?	
## output:
# t2: complete trajectory of teststatistic

strucchange.HL <- function(x,varmethod=c("window","acf","acfextra"),overlapping=TRUE,shiftcorrect=TRUE,borderN=10,...){
N <- length(x)

varmethod <- match.arg(varmethod)

	threedots <- list(...)
	index1 <- which(names(threedots) %in% c("dd","cc","K","type","order.max","aic","momentp"))
	if (length(index1)==0) threedots1 <- list() else{
		threedots1 <- threedots[index1]
		}
	index2 <- which(names(threedots) %in% c("type2","adjust","kernelused"))
	if (length(index2)==0) threedots2 <- list() else{
		threedots2 <- threedots[index2]
		}
	

if (varmethod=="window") {
	asy <- do.call(asymvar.window,c(list(x=x),list(overlapping=overlapping),list(shiftcorrect=shiftcorrect),list(borderN=borderN),list(obs="ranks"),threedots1))[[1]]
	}
if (varmethod=="acf") {
	asy <- do.call(asymvar.acf,c(list(x=x,shiftcorrect=shiftcorrect,obs="ranks"),threedots1))[[1]]	
	}
if (varmethod=="acfextra") {	
	asy <- do.call(asymvar.acfextra,c(list(x=x,shiftcorrect=shiftcorrect,obs="ranks"),threedots1))[[1]]
	}
t2 <- numeric(N-1)
for (m in 1:(N-1)){
	x2=x[(m+1):N]
	x1=x[1:m]
	n=N-m
	t1=meddiff(x1,x2)
		if (shiftcorrect) {
			x2 <- x2+t1
			}
		u0 <- do.call(densdiff,c(list(x1),list(x2),threedots2))
	t2[m]=N^(-3/2)*t1*m*n*u0	
	}
t2 <- t2/asy
return(t2)
}

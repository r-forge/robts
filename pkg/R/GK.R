###############
# auxiliary function: calculates (2-dim) correlation based on variances (GK-approach)
# input
# x: first variable to calculate correlation
# y: second variable to calculate correlation
# ...: parameters for Tau estimator
# method: what variance estimator to use (Qn, Tau and MAD available)
# output: correlation-coefficient
###############

GK <- function(x,y,method="Qn",...) {

### choose varaince estimator

if (!any(method==c("Qn","Tau","MAD"))) 
	{warning("This is no suitable variance estimator. Qn is used instead.")
	varest <- Qn}
if (method=="Qn") {
	varest <- Qn
	}
if (method=="Tau") {
	varest <- function(x) scaleTau2(x,...)
	}
if (method=="MAD") {
	varest <- mad
	}

### calculating necessary variances

varsum <- varest(x+y)^2
vardif <- varest(x-y)^2
if(varsum+vardif==0){warning("Something is wrong with your timeseries since variance estimation is 0.")
	return(NA)	
	}
return((varsum-vardif)/(varsum+vardif))
}
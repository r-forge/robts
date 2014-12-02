####################
# auxiliary function to make mediancor consistent
# input
# x: estimated mediancorrelation
# A: vector of simulated expection values of mediancorrelation under a true correlation of b
# b: vector of correlations were the expection value of mediancorrelation is simulated
####################

linearinterpol <- function(x,A,b) {

# the transformation to the sonsistent value is calculated by a linear interpolation based on the lattice of b respectivle A

ob <- which(x<=A)[1]	# determining the lower bound


# if the estimated correlation is not in the area where the inverse transformation is known


if(is.na(ob)){
	warning("The correlation is out of the simulated bounds. Therefore the correlation can not be transformed correctly.")
	return(max(b))
	}
if (ob==1){
	warning("The correlation is out of the simulated bounds. Therefore the correlation can not be transformed correctly.")
	return(min(b))
	}

un <- ob-1
wert <- b[un]+0.05/(A[ob]-A[un])*(x-A[un])	# linear interpolation
return(wert)
}
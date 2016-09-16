####################
# auxiliary function to make median correlation consistent
# input
# x: estimated median correlation
# a: vector of simulated expection values of mediancorrelation under a true correlation of b
# b: vector of correlations were the expection value of mediancorrelation is simulated
####################

linearinterpol <- function(x, a, b) {

  # the transformation to the sonsistent value is calculated by a linear interpolation based on the lattice of b respective a  
  ob <- which(x<=a)[1]	# determining the lower bound
  
  # check whether the estimated correlation is in the area where the inverse transformation is known:  
  if(is.na(ob)){
  	warning("The correlation is out of the simulated bounds. Therefore the correlation can not be transformed correctly.")
  	return(max(b))
  	}
  if (ob==1){
  	warning("The correlation is out of the simulated bounds. Therefore the correlation can not be transformed correctly.")
  	return(min(b))
  	}
  
  un <- ob - 1
  value <- b[un] + 0.05 / (a[ob] - a[un]) * (x - a[un])	# linear interpolation
  return(value)
}

###############
# auxiliary function: calculates bivariate correlation based on variances (GK-approach)
# input
# x: first variable to calculate correlation
# y: second variable to calculate correlation
# scalefn: what variance estimator to use (Qn, mad, scaleTau2, reweightedQN, ...)
# ...: arguments which are passed to the function scalefn
# output: correlation-coefficient
###############
  
  
corGK <- function(x, y, scalefn=Qn, ...) {
  #see the related function 'covGK' from package 'robustbase'
  
  # calculating necessary variances:  
  varsum <- scalefn(x+y, ...)^2
  vardif <- scalefn(x-y, ...)^2
  
  if(varsum+vardif==0){
    warning("Something is wrong since variance estimation is 0.")
  	return(NA)	
 	}
 	
 	result <- (varsum-vardif)/(varsum+vardif)
  return(result)
}

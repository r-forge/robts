####################
# auxiliary function: Huber-weight-function
# input
# x: vector of observations
# k: Robustness-Parameter (k=2.37 => 0.95 efficiency for locationestimation)
# output: robustified observation
#####################

Huber <- function(x,k=1.37) {
case <- abs(x) <= k 				# observation small enough to stay as it is?
return(x*case+k*sign(x)*(1-case))   # norm of larger values set to k 
}
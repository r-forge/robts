####################
# auxiliary function: Psi function
# input
# x: vector of observations
# k: Robustness-Parameter (k=1.37 => 0.95 efficiency for location estimation)
# output: robustified observation
#####################

M_psi <- function(x, type=c("huber", "bisquare"), k){
  type <- match.arg(type)
  if(missing(k)) k <- switch(type, huber=1.37, bisquare=4.68)
  if(type=="huber"){
    inner <- abs(x) <= k  # observation small enough to stay as it is?
    return(x*inner+k*sign(x)*(1-inner)) # norm of larger values set to k 
  } 
  if(type=="bisquare"){
    return(x*(1-(x/k)^2)^2*(abs(x)<=k))
  }
}

M_wgt <- function(x, type=c("huber", "bisquare"), k){
  type <- match.arg(type)
  if(missing(k)) switch(type, huber=1.37, bisquare=4.68)
  if(type=="huber"){
    return(apply(cbind(1, k/abs(x)), 1, min))
  } 
  if(type=="bisquare"){
    return(apply(cbind((3*k^4-3*x^2*k^2+x^4)/k^6,1/x^2), 1, min))
  }
}

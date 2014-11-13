####################
# auxiliary function: Tukeybiweight-function
# input
# x: vector of observations
# k: Robustness-Parameter (k=4.68 => 0.95 efficiency for locationestimation)
# output: robustified observation
#####################

Tukeyp <- function(x,k=4.68) {
return(x*(1-(x/k)^2)^2*(abs(x)<=k))   # usual tukey-weights
}
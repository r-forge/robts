######################################################################
##### Hodges-Lehmann Estimator of location difference ###############################
######################################################################

## calculates the median of all combinations of differences from first and second sample
## input:
# y: first sample
# z: second sample
# cor: additive correction (default is 0)
## output:
# estimator of difference between first and second sample

meddiff <- function(y, z, cor = 0){
  dif <- rep(y, each = length(z))
  dif <- dif - z
  med <- median(dif)
  res <- med - cor
  return(res)
}

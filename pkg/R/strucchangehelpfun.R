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

meddiff<-function(y,z,cor=0){
dif=rep(y,each=length(z))
dif=dif-z
med=median(dif)
return(med-cor)
}


## estimates the empirical distribution function at sample values
## input:
# y: sample
## output:
# erg: empirical distribution function at sample values
edf <- function(y) {
  T <- length(y)
  erg <- rank(y)/T
return(erg)
}



### quantileKS
## computes given Quantile of Kolmogorov-Smirnov distribution
# input: x in (0,1)
# output: q(x)

quantileKS <- function(x) {
x <- 1-x
pKS <- function(y) return(.Call(stats:::C_pKS2, p = y, 10^(-6))-x)
quan <- uniroot(pKS,lower=0.2,upper=3)
return(quan$root)
}


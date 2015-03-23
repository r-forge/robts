#### function calculates asymptotic confidence bounds for acf/pacf under white noise ####
## input:
# fun: what estimator is used
# n: sample size
# ci: covering probability
## output:
# lower an upper bound as vector

konfband <- function(fun,n,ci=0.95,...) {
quan <- qnorm(1-(1-ci)/2)
dotdotdot <- list(...)
if (fun=="acfGK") {
index <- which(names(dotdotdot)%in%"method")
if (length(index)==0) {
gr <- quan*sqrt(1/n*1/0.8227)
return(c(-gr,gr))
}
method <- dotdotdot[[index]]
if (method=="Qn") {
gr <- quan*sqrt(1/n*1/0.8227)
return(c(-gr,gr))
}
if (method=="Tau") {
gr <- quan*sqrt(1/n*1/0.8)
return(c(-gr,gr))
}
if (method=="MAD") {
gr <- quan*sqrt(1/n*1/0.3674)
return(c(-gr,gr))
}
}
if (fun=="acfrank") {
index <- which(names(dotdotdot)%in%"method")
if (length(index)==0) {
gr <- quan*sqrt(1/n*1)
return(c(-gr,gr))
}
if (method=="gaussian") {
gr <- quan*sqrt(1/n*1)
return(c(-gr,gr))
}
if (method=="spearman") {
gr <- quan*sqrt(1/n*pi^2/9)
return(c(-gr,gr))
}
if (method=="kendall") {
gr <- quan*sqrt(1/n*pi^2/9)
return(c(-gr,gr))
}
if (method=="quadrant") {
gr <- quan*sqrt(1/n*pi^2/4)
return(c(-gr,gr))
}

}
}




## estimates long run variance from acf
## input:
# y: timeseries
# cc: factor which determines, when autocorrelation is significant
# K: determines length of sequence in which at least one autocorrelation should be significantly large
# type: which kernal should be used
#	possibilities: 	"bartlett"
#			"trapezoid"
# tau: time of structural break if timeseries is not stationary
# h: jump height
## output:
# asy: estimated long run variance

asymvar.acf <- function(x, obs = c("untransformed", "ranks"), cc = 1.4, K = 3, type = c("bartlett", "trapezoid")){
  N <- length(x)
  type <- match.arg(type)
  obs <- match.arg(obs)
  if (obs=="ranks") {
    x <- rank(x)/length(x)
	}
  ac <- acf(x,plot=F,lag.max=floor(2*N**(2/3)))$acf
  vc <- cc*sqrt(log(N,base=10)/N)
  for (i in 1:floor(N**(2/3))){
    if(max(abs(ac[i+(1:K)]))<vc){break}
  }
  w <- ((2*i):1)/(2*i)
  if (type=="trapezoid"){w=c(rep(1,i),(i:1)/i)}
  if (sum(type==c("bartlett","trapezoid"))==0) {
    warning(paste(type,"is not implemented. Using bbartlett kernal instead."))
  }
  ac <- acf(x,plot=F,type="covariance",lag.max=2*i)$acf[1:(2*i+1)]
  asy <- ac[1]+2*sum(ac[2:(2*i+1)]*w)
  erg <- list(lrv=asy,bandwith=i)
  return(erg)
}

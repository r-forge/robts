arrob <- function(x, aic = TRUE, order.max, method = c("yw", "regression", "gm", "filter"), na.action = na.fail, series = deparse(substitute(x)),asyvar=FALSE,bootpar=list(num=100,blockl=floor(2*length(x)^(1/2))), ...) {

  # checks and preparations:
	method <- match.arg(method)
	# more checks are done by the function of the respective method

	res <- switch(method,
	 yw = arrob.yw(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 regression = arrob.regression(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 gm = arrob.gm(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...),
	 filter = arrob.filter(x=x, aic=aic, order.max=order.max, na.action=na.action, series=series, ...)
  )
  res$call <- match.call()
if (asyvar==TRUE)	{
if (res$order==0) return(res)
n <- length(x)
blockl <- bootpar$blockl
num <- bootpar$num
m <- ceiling(n/blockl)
startv <- matrix(sample(1:(n-blockl+1),replace=TRUE,m*num),ncol=num)
tsmat <- apply(startv,2,fuwo,x=x,l=blockl)[1:n,]

if (method=="yw") fumu <- function(xx) return(arrob.yw(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)
if (method=="regression") fumu <- function(xx) return(arrob.regression(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)
if (method=="gm") fumu <- function(xx) return(arrob.gm(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)
if (method=="filter") fumu <- function(xx) return(arrob.filter(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)

bootest <- t(apply(tsmat,2,fumu))
if(res$order==1) res$asy.var.coef <- var(t(bootest)) else res$asy.var.coef <- cov(bootest)

}
  
	return(res)
}

fuwo <- function(x,st,l) {
y <- numeric(length(st)*l)
for(i in 1:length(st)) y[((i-1)*l+1):(i*l)] <- x[st[i]:(st[i]+l-1)] 
return(y)
}

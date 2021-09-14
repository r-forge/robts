armarob <- function(x, arorder, maorder, aic=FALSE, aicpenalty = function(p) {2*p}, na.action = na.fail, series = deparse(substitute(x)),asyvar=FALSE,bootpar=list(num=100,blockl=floor(2*length(x)^(1/2))), ...) {
n <- length(x)
xfreq <- frequency(x)
arest <- arrob.filter(x,order.max=4*(arorder+maorder),aic=TRUE)
if (aic==FALSE) {
    arorderv <- arorder
    maorderv <- maorder
    if (length(arest$ar)==0) {startpar <- c(rep(0,arorder+maorder))} else{

        if (arorder>0) {
            acf1 <- ARMAacf(ar=arest$ar,lag.max=arorder)[(1:arorder)+1]
            acf2 <- numeric(maorder)
            for (i in 1:maorder) {
                acf2[i] <- cov(arest$x[(i+1):n],M_psi(arest$resid[1:(n-i)]/sqrt(arest$var.pred), type = "smooth", k = c(2, 3)),use="complete.obs")/sqrt(arest$var.pred)*1.184
                }
            fu <- function(par) {
                res <- numeric(arorder+maorder)
                for (i in 1:arorder) {
                    res[i] <- ARMAacf(ar=par[1:arorder],ma=par[(arorder+1):(arorder+maorder)],lag.max=i)[i+1]-acf1[i]
                    }
                for (i in 1:maorder) {
                    res[arorder+i] <- acfu(i,ar=par[1:arorder],ma=par[(arorder+1):(arorder+maorder)])-acf2[i]
                    }
                return(sum(abs(res)))
                }
            startpar <- optim(rep(0,arorder+maorder),fu)$par
        } else {

            acf2 <- numeric(maorder)
            for (i in 1:maorder) {
                acf2[i] <- cov(arest$x[(i+1):n],M_psi(arest$resid[1:(n-i)]/sqrt(arest$var.pred), type = "smooth", k = c(2, 3)),use="complete.obs")/sqrt(arest$var.pred)*1.184
                }
            startpar <- acf2
            }

        }
    if(arorder>0) {
        fu2 <- function(par) {
            residualv <- filterrob.statespaceARMA(x, par[1:arorder], par[(arorder+1):(arorder+maorder)], arest$var.pred, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals[(1+max(arorder,maorder+1)):n]
            return(scaleTau2(residualv))
            }

        estpar <- optim(startpar,fu2)
        filtered <- filterrob.statespaceARMA(x, estpar$par[1:arorder], estpar$par[(arorder+1):(arorder+maorder)], estpar$value, median, psi.l = 2, psi.0 = 3, na.action = na.fail)
        RAICS <- log(estpar$value)+ aicpenalty(arorder+maorder+1)/(n-arorder-maorder)
        var.predv <- estpar$value
        arestv <- estpar$par[1:arorder]
        maestv <- estpar$par[(arorder+1):(arorder+maorder)]
        residualv <- filtered$residuals
        } else{
            fu2 <- function(par) {
                residualv <- filterrob.statespaceARMA(x, 0, par, arest$var.pred, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals[(maorder+2):n]
                return(scaleTau2(residualv))
                }
            if (length(startpar)>1) {
                estpar <- optim(startpar,fu2)
                filtered <- filterrob.statespaceARMA(x, 0, estpar$par, estpar$value, median, psi.l = 2, psi.0 = 3, na.action = na.fail)
                RAICS <- log(estpar$value)+ aicpenalty(maorder+1)/(n-maorder)
                var.predv <- estpar$value
                arestv <- NULL
                maestv <- estpar$par
                residualv <- filtered$residuals
            } else {
                estpar <- nlm(fu2,startpar)  
                filtered <- filterrob.statespaceARMA(x, 0, estpar$estimate, estpar$minimum, median, psi.l = 2, psi.0 = 3, na.action = na.fail)
                RAICS <- log(estpar$minimum)+ aicpenalty(maorder+1)/(n-maorder) 
                var.predv <- estpar$minimum
                arestv <- NULL
                maestv <- estpar$estimate 
                residualv <- filtered$residuals                  
                }
            
            }

    } else {
        aictable <- matrix(ncol=maorder+1,nrow=arorder+1)
        aictable[,1] <- arest$aic[1:(arorder+1)]
        for (i1 in 1:maorder) {
            for (j1 in 0:arorder) {
              
                if (j1>0) {
                    acf1 <- ARMAacf(ar=arest$ar,lag.max=j1)[(1:j1)+1]
                    acf2 <- numeric(i1)
                    for (i in 1:i1) {
                        acf2[i] <- cov(arest$x[(i+1):n],M_psi(arest$resid[1:(n-i)]/sqrt(arest$var.pred), type = "smooth", k = c(2, 3)),use="complete.obs")/sqrt(arest$var.pred)*1.184
                        }
                    fu <- function(par) {
                        res <- numeric(j1+i1)
                        for (i in 1:j1) {
                            res[i] <- ARMAacf(ar=par[1:j1],ma=par[(j1+1):(j1+i1)],lag.max=i)[i+1]-acf1[i]
                            }
                        for (i in 1:i1) {
                            res[j1+i] <- acfu(i,ar=par[1:j1],ma=par[(j1+1):(j1+i1)])-acf2[i]
                            }
                    return(sum(abs(res)))
                    }
                startpar <- optim(rep(0,arorder+maorder),fu)$par
                } else {

                    acf2 <- numeric(i1)
                    for (i in 1:i1) {
                        acf2[i] <- cov(arest$x[(i+1):n],M_psi(arest$resid[1:(n-i)]/sqrt(arest$var.pred), type = "smooth", k = c(2, 3)),use="complete.obs")/sqrt(arest$var.pred)*1.184
                        }
                    startpar <- acf2
                    }
            
                if(j1>0) {
                    fu2 <- function(par) {
                        residualv <- filterrob.statespaceARMA(x, par[1:j1], par[(j1+1):(j1+i1)], arest$var.pred, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals[(1+max(j1,i1+1)):n]
                        return(scaleTau2(residualv))
                        }

                    estpar <- optim(startpar,fu2)
                    aicv <- log(estpar$value)+ aicpenalty(j1+i1+1)/(n-j1-i1)
                    if (min(aictable,na.rm=TRUE) > aicv) 
                        {
                        arestv <- estpar$par[1:j1]
                        maestv <- estpar$par[(j1+1):(j1+i1)]
                        var.predv <- estpar$value

                        }
                    aictable[j1+1,i1+1] <- aicv
                    } else{
                    fu2 <- function(par) {
                        residualv <- filterrob.statespaceARMA(x, 0, par, arest$var.pred, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals[(i1+2):n]
                        return(scaleTau2(residualv))
                        }
                    if (length(startpar)>1) {
                        estpar <- optim(startpar,fu2)
                        aicv <- log(estpar$value)+ aicpenalty(i1+1)/(n-i1)
                        if (min(aictable,na.rm=TRUE) > aicv) 
                            {
                            arestv <- NULL
                            maestv <- estpar$par[1:i1]
                            var.predv <- estpar$value


                            }
                        aictable[j1+1,i1+1] <- aicv
                        
                        } else {
                        estpar <- nlm(fu2,startpar)
                        aicv <- log(estpar$minimum)+ aicpenalty(i1+1)/(n-i1)
                        if (min(aictable,na.rm=TRUE) > aicv) 
                            {
                            arestv <- NULL
                            maestv <- estpar$estimate[1:i1]
                            var.predv <- estpar$minimum

    
                            }
                        aictable[j1+1,i1+1] <- aicv
                     
                        }
                  
                    }
            }
        }
    rownames(aictable) <- paste("AR-order ", 1:3,":",sep="")
    colnames(aictable) <- paste("MA-order ", 1:3,":",sep="")
    if ((arorderv==0)&(maorderv>0)) residualv <- filterrob.statespaceARMA(x, 0, maestv, var.predv, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals
    if ((arorderv>0)&(maorderv>0)) residualv <- filterrob.statespaceARMA(x, arestv, maestv, var.predv, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals
    if ((arorderv=0)&(maorderv=0)) residualv <- x-median(x)
    if ((arorderv>0)&(maorderv==0)) residualv <- filterrob.statespace(x, arestv, var.predv, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals
    }

  res <- list(
		order = c(length(arestv),length(maestv)),
		ar = arestv,
        ma = maestv,
		var.pred = var.predv,
		x.mean = median(filtered$filtered),
		x.intercept = NULL,
        aic = aictable
		n.used = n,
		resid = residualv,
		method = "filter",
		series = series,
		frequency = xfreq,
		asy.var.coef = NULL,
		x = x
	)
	attr(res, "na.action") <- attr(x, "na.action")
	class(res) <- c("arrob", "arma")
	return(res)
}
#if (asyvar==TRUE)	{
#if (res$order==0) return(res)
#n <- length(x)
#blockl <- bootpar$blockl
#num <- bootpar$num
#m <- ceiling(n/blockl)
#startv <- matrix(sample(1:(n-blockl+1),replace=TRUE,m*num),ncol=num)
#tsmat <- apply(startv,2,fuwo,x=x,l=blockl)[1:n,]

#if (method=="yw") fumu <- function(xx) return(arrob.yw(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)
#if (method=="regression") fumu <- function(xx) return(arrob.regression(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)
#if (method=="mm") fumu <- function(xx) return(arrob.gm(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)
#if (method=="filter") fumu <- function(xx) return(arrob.filter(xx,aic=FALSE,order.max=res$order,na.action=na.action,...)$ar)

#bootest <- t(apply(tsmat,2,fumu))
#if(res$order==1) res$asy.var.coef <- var(t(bootest)) else res$asy.var.coef <- cov(bootest)

#}
  
#	return(res)
#}

#fuwo <- function(x,st,l) {
#y <- numeric(length(st)*l)
#for(i in 1:length(st)) y[((i-1)*l+1):(i*l)] <- x[st[i]:(st[i]+l-1)] 
#return(y)
#}

acfu <- function(d,ar,ma) {
p <- length(ar)
q <- length(ma)
if (d > p) ar <- c(ar,rep(0,d-p))
if (d > q) ma <- c(ma,rep(0,d-q))
v <- 0
arv <- rep(0,d+1)
arv[1] <- 1
for (i in 1:d) {
arv[(1+i):(d+1)] <- arv[i]*ar[1:(d-i+1)]+arv[(1+i):(d+1)]
}
for (i in 1:d) {
v <- v+arv[i]*ma[d-i+1]
}
v <- v+arv[d+1]
return(v)
} 

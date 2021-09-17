armarob <- function(x, arorder, maorder, aic=FALSE, aicpenalty = function(p) {2*p}, na.action = na.fail, series = deparse(substitute(x))) {
n <- length(x)
xfreq <- frequency(x)
arest <- arrob.filter(x,order.max=4*(arorder+maorder),aic=TRUE)
if (aic==FALSE) {
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
                    tryacf <- try(tacvfARMA(phi=par[1:arorder],theta=-par[(arorder+1):(arorder+maorder)],maxLag=i),silent=TRUE)
                    if (inherits(tryacf,"try-error")) return(Inf)
                    res[i] <- tryacf[i+1]/tryacf[1]-acf1[i]
                    }
                for (i in 1:maorder) {
                    res[arorder+i] <- acfu(i,ar=par[1:arorder],ma=par[(arorder+1):(arorder+maorder)])-acf2[i]
                    }
                return(sum(abs(res)))
                }
            arstart <- rep(0,arorder)
            if (arest$order>0) arstart[1:min(arest$order,arorder)] <- arest$ar[1:min(arest$order,arorder)]
            startpar <- optim(c(arstart,rep(0,maorder)),fu)$par
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
            residualv <- try(filterrob.statespaceARMA(x, par[1:arorder], par[(arorder+1):(arorder+maorder)], arest$var.pred, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals[(1+max(arorder,maorder+1)):n],silent=TRUE)
            if(inherits(residualv,"try-error")) return(Inf)
            return(scaleTau2(residualv))
            }

        estpar <- try(optim(startpar,fu2),silent=TRUE)
        if(inherits(estpar,"try-error")) {estpar <- optim(rep(0,arorder+maorder),fu2)}
        filtered <- filterrob.statespaceARMA(x, estpar$par[1:arorder], estpar$par[(arorder+1):(arorder+maorder)], estpar$value, median, psi.l = 2, psi.0 = 3, na.action = na.fail)
        aictable <- 2*log(estpar$value)+ aicpenalty(arorder+maorder+1)/(n-arorder-maorder)
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
                aictable <- 2*log(estpar$value)+ aicpenalty(maorder+1)/(n-maorder)
                var.predv <- estpar$value
                arestv <- NULL
                maestv <- estpar$par
                residualv <- filtered$residuals
            } else {
                estpar <- nlm(fu2,startpar)  
                filtered <- filterrob.statespaceARMA(x, 0, estpar$estimate, estpar$minimum, median, psi.l = 2, psi.0 = 3, na.action = na.fail)
                aictable <- 2*log(estpar$minimum)+ aicpenalty(maorder+1)/(n-maorder) 
                var.predv <- estpar$minimum
                arestv <- NULL
                maestv <- estpar$estimate 
                residualv <- filtered$residuals                  
                }
            
            }

    } else {
        aictable <- matrix(ncol=maorder+1,nrow=arorder+1)
        aictable[,1] <- arest$aic[1:(arorder+1)]
        maestv <- NULL
        arestv <- arest$ar
        var.predv <- arest$var.pred
        
        for (i1 in 1:maorder) {
            for (j1 in 0:arorder) {
              
                if (j1>0) {
                    if(length(arest$ar)>0) {
                        acf1 <- ARMAacf(ar=arest$ar,lag.max=j1)[(1:j1)+1]
                        } else {acf1 <- rep(0,j1)} 
    
                    acf2 <- numeric(i1)
                    for (i in 1:i1) {
                        acf2[i] <- cov(arest$x[(i+1):n],M_psi(arest$resid[1:(n-i)]/sqrt(arest$var.pred), type = "smooth", k = c(2, 3)),use="complete.obs")/sqrt(arest$var.pred)*1.184
                        }
                    fu <- function(par) {
                        res <- numeric(j1+i1)
                        for (i in 1:j1) {
                            tryacf <- try(tacvfARMA(phi=par[1:j1],theta=-par[(j1+1):(j1+i1)],maxLag=i),silent=TRUE)
                            if (inherits(tryacf,"try-error")) return(Inf)
                            res[i] <- tryacf[i+1]/tryacf[1]-acf1[i]
                            }
                        for (i in 1:i1) {
                            res[j1+i] <- acfu(i,ar=par[1:j1],ma=par[(j1+1):(j1+i1)])-acf2[i]
                            }
                    return(sum(abs(res)))
                    }
                arstart <- rep(0,j1)
                mastart <- rep(0,i1)
                if (length(arestv)>0) arstart[1:min(j1,length(arestv))] <- arestv[1:min(j1,length(arestv))]
                if (length(maestv)>0) mastart[1:min(j1,length(maestv))] <- maestv[1:min(j1,length(maestv))]
                startpar <- optim(c(arstart,mastart),fu)$par
                } else {

                    acf2 <- numeric(i1)
                    for (i in 1:i1) {
                        acf2[i] <- cov(arest$x[(i+1):n],M_psi(arest$resid[1:(n-i)]/sqrt(arest$var.pred), type = "smooth", k = c(2, 3)),use="complete.obs")/sqrt(arest$var.pred)*1.184
                        }
                    startpar <- acf2
                    }
            
                if(j1>0) {
                    fu2 <- function(par) {
                        residualv <- try(filterrob.statespaceARMA(x, par[1:j1], par[(j1+1):(j1+i1)], arest$var.pred, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals[(1+max(j1,i1+1)):n],silent=TRUE)
                        if (inherits(residualv,"try-error")) return(Inf)
                        return(scaleTau2(residualv))
                        }

                    try(estpar <- optim(startpar,fu2),silent=TRUE)
                    if (inherits(estpar,"try-error")) {aicv <- Inf} else {
                    aicv <- 2*log(estpar$value)+ aicpenalty(j1+i1+1)/(n-j1-i1)}
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
                        aicv <- 2*log(estpar$value)+ aicpenalty(i1+1)/(n-i1)
                        if (min(aictable,na.rm=TRUE) > aicv) 
                            {
                            arestv <- NULL
                            maestv <- estpar$par[1:i1]
                            var.predv <- estpar$value


                            }
                        aictable[j1+1,i1+1] <- aicv
                        
                        } else {
                        estpar <- nlm(fu2,startpar)
                        aicv <- 2*log(estpar$minimum)+ aicpenalty(i1+1)/(n-i1)
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
    rownames(aictable) <- paste("AR-order ", 0:2,":",sep="")
    colnames(aictable) <- paste("MA-order ", 0:2,":",sep="")
    if ((length(arestv)==0)&(length(maestv)>0)) residualv <- filterrob.statespaceARMA(x, 0, maestv, var.predv, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals
    if ((length(arestv)>0)&(length(maestv)>0)) residualv <- filterrob.statespaceARMA(x, arestv, maestv, var.predv, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals
    if ((length(arestv)==0)&(length(maestv)==0)) residualv <- x-median(x)
    if ((length(arestv)>0)&(length(maestv)==0)) residualv <- filterrob.statespace(x, arestv, var.predv, median, psi.l = 2, psi.0 = 3, na.action = na.fail)$residuals
    }

  res <- list(
		order = c(length(arestv),length(maestv)),
		ar = arestv,
        ma = maestv,
		var.pred = var.predv^2,
		x.mean = median(x),
		x.intercept = NULL,
        aic = aictable,
		n.used = n,
		resid = residualv,
		method = "filter",
		series = series,
		frequency = xfreq,
		asy.var.coef = NULL,
		x = x
	)
	attr(res, "na.action") <- attr(x, "na.action")
	class(res) <- c("armarob")
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

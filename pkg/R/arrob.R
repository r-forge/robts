arrob <- function(x, aic = TRUE, order.max = NULL,
	method = c("yule-walker", "durbin-levinson", "robustregression", "filter", "gm"),
	na.action = na.fail, series = deparse(substitute(x)), ...,
	acf.approach = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim", "acfrobfil", "acfregression"), aicpenalty=function(p) {2*p}) {
	cl <- match.call()
	method <- match.arg(method)
	x <- na.action(as.ts(x))
	
	acf.approach = match.arg(acf.approach)
	
	x <- na.action(as.ts(x))
	n <- length(x)
	if (is.null(order.max)) order.max <- min(c((n - 1) / 4, 10 * log(n, 10)))
	if (aic) {
		if (method == "yule-walker") {
			re <- ARopt.YW(x, pmax = order.max, acf.approach = acf.approach,aicpenalty=aicpenalty,...)
			if (sum(is.na(re))>0) {
				warning("Calculation failed")
				return(NA)
			}
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			residuals <- re$residuals
			x.mean <- re$x.mean
		}
		if (method == "durbin-levinson") {
			re <- ARopt.acf(tss = x, pmax = order.max, acf.approach = acf.approach,aicpenalty=aicpenalty,...)
			if (sum(is.na(re))>0) {
				warning("Calculation failed")
				return(NA)
			}
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			residuals <- re$residuals
			x.mean <- re$x.mean
		}
		if (method == "robustregression") {
			re <- lmrobARopt(x, pmax = order.max, interc = TRUE,aicpenalty=aicpenalty, ...)
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			residuals <- re$resid
			x.mean <- re$x.mean
		}
		if (method == "filter") {
			re <- ARopt.filter(x, pmax = order.max,aicpenalty=aicpenalty, ...)
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			residuals <- re$resid
			x.mean <- re$x.mean
		}
		if (method == "gm") {
			re <- bestAR(x, maxp = order.max,aicpenalty=aicpenalty, ...)
			wm <- which.min(re$aic)[1]
			if (wm == 1) ph <- NULL else ph <- re$phimatrix[wm, 1:(wm - 1)]
			x.mean=re$x.mean
			aicv <- re$aic
			var.pred <- re$var.pred[wm]
			residuals <- re$residuals[,wm]
			if (wm != 1) residuals <- residuals[-(1:(wm - 1))]
			partialacf <- diag(re$phimatrix[-1,])
		}
	} else {
		if (method == "yule-walker") {
			acorf <- as.numeric(acfrob(x, approach = acf.approach, lag.max=order.max,plot = FALSE, ...)$acf)
			ph <- solveYW(acorf, p = order.max)
			x.mean <- median(x)
			D <- matrix(nrow = n - order.max, ncol = order.max)
			for (j in 1:order.max) D[, j] <- x[(order.max + 1 - j):(n - j)]
			D <- cbind(1, D)
			ph1 <- c(x.mean * (1 - sum(ph)), ph)
			residuals <- x[(order.max + 1):n] - D %*% ph1
			var.pred <- Qn(residuals)^2
			partialacf <- ARMAacf(ar=ph,lag.max=order.max,pacf=TRUE)
		}
		if (method == "durbin-levinson") {
			acorf <- as.numeric(acfrob(x, approach = acf.approach, lag.max=order.max,plot = FALSE, ...)$acf)
			ph <- DurbinLev(acorf)[[1]][[order.max]]
			x.mean <- median(x)
			D <- matrix(nrow = n - order.max, ncol = order.max)
			for (j in 1:order.max) D[, j] <- x[(order.max + 1 - j):(n - j)]
			D <- cbind(1, D)
			ph1 <- c(x.mean * (1 - sum(ph)), ph)
			residuals <- x[(order.max + 1):n] - D %*% ph1
			var.pred <- Qn(residuals)^2
			partialacf <- ARMAacf(ar=ph,lag.max=order.max,pacf=TRUE)
		}
		if (method == "robustregression") {
			re <- lmrobAR(x, p = order.max, interc = TRUE, ...)
			x.mean <- re$coefficients[order.max+1]
			ph <- re$coefficients[1:order.max]
			var.pred <- re$model$scale^2
			residuals <- re$model$residuals
			partialacf <- ARMAacf(ar=ph,lag.max=order.max,pacf=TRUE)
		}
		if (method == "filter") {
			re <- ARfilter(x, p = order.max, ...)
			ph <- re[[6]][order.max,]
			x.mean <- scaleTau2(x,mu.too=TRUE)[1]
			var.pred <- re[[2]][order.max]^2
			partialacf <- re[[1]]
			xcen <- re[[5]][,order.max]-x.mean
			xcenma <- matrix(ncol=order.max,nrow=n-order.max)
			for (i in 1:order.max) {xcenma[,i] <- xcen[(order.max-i+1):(n-i)]}
			residuals <- x[(order.max+1):n]-x.mean-xcenma%*%ph
		}
		if (method == "gm") {
			re <- bestAR(x, maxp = order.max, ...)
			ph <- re$phimatrix[order.max+1, 1:order.max]
			x.mean=re$x.mean
			var.pred <- re$var.pred[order.max+1]
			residuals <- re$residuals[-(1:order.max), order.max + 1]
			partialacf <- diag(re$phimatrix[-1,])
		}
	}
	p <- length(ph)
	if(!aic) aicv <- log(var.pred) + aicpenalty(order.max+1)/(n-order.max)
	class(residuals) <- class(x)
	res <- list(
		order = p,
		ar = ph,
		var.pred = var.pred,
		x.mean = x.mean,
		x.intercept = NULL,
		aic = aicv,
		n.used = n,
		order.max = order.max,
		partialacf = partialacf,
		resid = residuals,
		method = method,
		series = series,
		frequency = frequency(x),
		call = cl,
		asy.var.coef = NULL
	)
	class(res) <- "ar"
	return(res)
}

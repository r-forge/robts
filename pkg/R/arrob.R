arrob <- function(x, aic = TRUE, order.max = NULL,
	method = c("yule-walker", "durbin-levinson", "robustregression", "filter", "gm"),
	na.action = na.fail, series = deparse(substitute(x)), ...,
	acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"),aicpenalty=function(p) {return(2*p)}) {
	
	method <- match.arg(method)
	if (!any(method == c("yule-walker", "durbin-levinson", "robustregression", "filter", "gm"))) stop("No valid method chosen.")
	acf.fun = match.arg(acf.fun)
	if (!any(acf.fun == c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"))) stop("No valid acf function chosen.")
	
	x <- na.action(x)
	n <- length(x)
	if (is.null(order.max)) order.max <- min((n - 1) / 4, 10 * log(n, 10))
	
	if (aic) {
		if (method == "yule-walker") {
			re <- ARopt.YW(x, pmax = order.max, acf.fun = acf.fun,aicpenalty=aicpenalty)
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			resid <- re$resid
			x.mean <- re$x.mean
		}
		if (method == "durbin-levinson") {
			re <- ARopt.acf(tss = x, pmax = order.max, acf.fun = acf.fun,aicpenalty=aicpenalty)
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			resid <- re$resid
			x.mean <- re$x.mean
		}
		if (method == "robustregression") {
			re <- lmrobARopt(x, pmax = order.max, interc = TRUE,aicpenalty=aicpenalty, ...)
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			resid <- re$resid
			x.mean <- re$x.mean
		}
		if (method == "filter") {
			re <- ARopt.filter(x, pmax = order.max,aicpenalty=aicpenalty, ...)
			aicv <- re$aic
			ph <- re$coefficients
			var.pred <- re$var.pred
			partialacf <- re$partialacf
			resid <- re$resid
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
			partialacf <- diag(re$phimatrix[-1,])
		}
	} else {
		if (method == "yule-walker") {
			acorf <- as.numeric(acfrob(x, fun = acf.fun, lag.max=order.max,plot = FALSE, ...)$acf)
			ph <- solveYW(acorf, p = order.max)
			x.mean <- median(x)
			D <- matrix(nrow = n - order.max, ncol = order.max)
			for (j in 1:order.max) D[, j] <- x[(order.max + 1 - j):(n - j)]
			D <- cbind(1, D)
			ph1 <- c(x.mean * (1 - sum(ph)), ph)
			residuals <- x[(order.max + 1):n] - D %*% ph1
			var.pred <- Qn(residuals)^2
			residuals <- c(rep(NA,order.max),residuals)
			partialacf <- ARMAacf(ar=ph,lag.max=order.max,pacf=TRUE)
		}
		if (method == "durbin-levinson") {
			acorf <- as.numeric(acfrob(x, fun = acf.fun, lag.max=order.max,plot = FALSE, ...)$acf)
			ph <- DurbinLev(acorf)[[1]][[order.max]]
			x.mean <- median(x)
			D <- matrix(nrow = n - order.max, ncol = order.max)
			for (j in 1:order.max) D[, j] <- x[(order.max + 1 - j):(n - j)]
			D <- cbind(1, D)
			ph1 <- c(x.mean * (1 - sum(ph)), ph)
			residuals <- x[(order.max + 1):n] - D %*% ph1
			var.pred <- Qn(residuals)^2
			residuals <- c(rep(NA,order.max),residuals)
			partialacf <- ARMAacf(ar=ph,lag.max=order.max,pacf=TRUE)
		}
		if (method == "robustregression") {
			re <- lmrobAR(x, p = order.max, interc = TRUE, ...)
			x.mean <- re$coefficients[order.max+1]
			ph <- re$coefficients[1:order.max]
			var.pred <- re$model$scale
			residuals <- c(rep(NA,order.max),re$model$residuals)
			partialacf <- ARMAacf(ar=ph,lag.max=order.max,pacf=TRUE)
		}
		if (method == "filter") {
			re <- ARfilter(x, p = order.max, ...)
			ph <- re[[6]][order.max,]
			x.mean <- scaleTau2(x,mu.too=TRUE)[1]
			var.pred <- re[[2]][order.max]
			partialacf <- re[[1]]
			xcen <- re[[5]][,order.max]-x.mean
			xcenma <- matrix(ncol=order.max,nrow=n-order.max)
			for (i in 1:order.max) {xcenma[,i] <- xcen[(order.max-i+1):(n-i)]}
			residuals <- c(rep(NA,order.max),x[(order.max+1):n]-x.mean-xcenma%*%ph)
		}
		if (method == "gm") {
			re <- bestAR(x, maxp = order.max, ...)
			ph <- re$phimatrix[order.max+1, 1:order.max]
			x.mean=re$x.mean
			var.pred <- re$var.pred[order.max+1]
			residuals <- re$residuals[,order.max+1]
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
		call = NULL,
		asy.var.coef = NULL
	)
	class(res) <- "ar"
	return(res)
}

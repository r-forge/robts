arrob <- function(x, aic = TRUE, order.max = NULL,
	method = c("yule-walker", "durbin-levinson", "ols", "filter"),
	na.action = na.fail, series = deparse(substitute(x)), ...,
	acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim")) {
	
	method <- match.arg(method)
	if (!any(method == c("yule-walker", "durbin-levinson", "ols", "filter"))) stop("No valid method chosen.")
	acf.fun = match.arg(acf.fun)
	if (!any(acf.fun == c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"))) stop("No valid acf function chosen.")
	
	x <- na.action(x)
	n <- length(x)
	if (is.null(order.max)) order.max <- min((n - 1) / 4, 10 * log(n, 10))
	
	if (aic) {
		if (method == "yule-walker") {
			ph <- ARopt.YW(x, pmax = order.max, acf.fun = acf.fun)
		}
		if (method == "durbin-levinson") {
			ph <- ARopt.acf(tss = x, aic = TRUE, pmax = order.max, acf.fun = acf.fun)
		}
		if (method == "ols") {
			ph <- lmrobARopt(x, pmax = order.max, interc = FALSE, ...)$coefficients
		}
		if (method == "filter") {
			popt <- ARopt.filter(x, pmax = order.max, ...)
			acorf <- ARfilter(x, p = popt, ...)[[5]]
			ph <- solveYW(acorf, p = popt)
		}
	} else {
		if (method == "yule-walker") {
			acorf <- as.numeric(acfrob(x, fun = acf.fun, plot = FALSE, ...)$acf)
			ph <- solveYW(acorf, p = order.max)
		}
		if (method == "durbin-levinson") {
			ph <- ARopt.acf(tss = x, aic = FALSE, pmax = order.max, acf.fun = acf.fun)
		}
		if (method == "ols") {
			ph <- lmrobAR(x, p = order.max, interc = FALSE, ...)$coefficients
		}
		if (method == "filter") {
			acorf <- ARfilter(x, p = order.max, ...)[[5]]
			ph <- solveYW(acorf, p = order.max)
		}
	}
	p <- length(ph)
	# residuals:
	D <- matrix(nrow = n - p, ncol = p)
	for (j in 1:p) D[, j] <- x[(p + 1 - j):(n - j)]
	D <- cbind(1, D)
	ph1 <- c(median(x) * (1 - sum(ph)), ph)
	resid <- x[(p + 1):n] - D %*% ph1
	resid <- c(rep(NA, p), resid)
	class(resid) <- class(x)
	res <- list(
		order = p,
		ar = ph,
		var.pred = NULL,
		x.mean = median(x),
		x.intercept = NULL,
		aic = NULL,
		n.used = n,
		order.max = order.max,
		partialacf = NULL,
		resid = resid,
		method = method,
		series = series,
		frequency = frequency(x),
		call = NULL,
		asy.var.coef = NULL
	)
	class(res) <- "ar"
	return(res)
}
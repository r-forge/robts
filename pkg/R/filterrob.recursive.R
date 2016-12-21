filterrob.recursive <- function(x, ar, var.pred, locfn, psi.l = 2, psi.0 = 3, na.action = na.fail) {
  ists <- is.ts(x)
  if (ists) xtsp <- tsp(x)
  x <- handle_missings_ts(x, na.action)
  if (!is.numeric(x)) stop("'x' must be numeric")
  psifn <- function(x) M_psi(x, type = "smooth", k = c(psi.l, psi.0))
  if(missing(locfn)) locfn <- function(x) scaleTau2(x, mu.too = TRUE)[1]
  p <- length(ar)
  mu <- locfn(x)
  sigma <- sqrt(var.pred)
  xfilt <- c(rep.int(0, p), x - mu)
  xhat <- numeric(length(x))
  residuals <- numeric(length(x))
  for (i in seq_along(x)) {
    xhat[i] <- ifelse(p > 0, sum(ar * xfilt[p + i - 1:p]), 0)
    residuals[i] <- xfilt[p + i] - xhat[i]
    xfilt[p + i] <- xhat[i] + psifn(residuals[i]/sigma) * sigma
  }
  residuals[seq_len(p)] <- NA
  filtered <- xfilt[p + seq_along(x)] + mu
  residuals <- naresid(attr(x, "na.action"), residuals)
  filtered <- napredict(attr(x, "na.action"), filtered)  
  if (ists) {
    attr(residuals, "tsp") <- attr(filtered, "tsp") <- xtsp
    attr(residuals, "class") <- attr(filtered, "class") <- "ts"
  } 
  res <- list(filtered = filtered, residuals = residuals)
 	attr(res, "na.action") <- attr(x, "na.action")
  return(res)
}

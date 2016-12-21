residuals.arrob <- function(object, x = object$x, method = c("recursive", "statespace", "nonrobust"), na.action = na.fail, ...) {
  method <- match.arg(method)
  if(method == "nonrobust"){ #nonrobust filtering 
    filterout <- filterrob(x = x, ar = object$ar, var.pred = object$var.pred, method = "recursive", psi.l = 1e+3, psi.0 = 1e+4, na.action = na.action)
  } else { 
    filterout <- filterrob(x = x, ar = object$ar, var.pred = object$var.pred, method = method, ..., na.action = na.action)
  }
  res <- naresid(attr(x, "na.action"), filterout$residuals)
  attr(res, "na.action") <- attr(filterout, "na.action")
  return(res)
}

fitted.arrob <- function(object, x = object$x, method = c("recursive", "statespace", "nonrobust"), na.action = na.fail, ...) {
  residuals <- residuals.arrob(object = object, x = x, method = method, na.action = na.action, ...)
  res <- object$x - residuals
  res <- napredict(attr(x, "na.action"), res)
  attr(res, "na.action") <- attr(residuals, "na.action") 
  return(res)
}

filtered <- function(object, ...) UseMethod("filtered")

filtered.arrob <- function(object, x = object$x, method = c("recursive", "statespace", "nonrobust"), na.action, ...) {
  method <- match.arg(method)
  if(method == "nonrobust"){  
    filterout <- filterrob(x = x, ar = object$ar, var.pred = object$var.pred, method = "recursive", psi.l = 1e+3, psi.0 = 1e+4, na.action = na.action) # nonrobust filtering, currently realized by chosing the tuning constants to be very large
  } else { 
    filterout <- filterrob(x = x, ar = object$ar, var.pred = object$var.pred, method = method, na.action = na.action, ...)
  }
  res <- napredict(attr(x, "na.action"), filterout$filtered)
  attr(res, "na.action") <- attr(filterout, "na.action")
  return(res)
}

predict.arrob <- function(object, newdata = object$x, n.ahead = 1, se.fit = TRUE, method = c("recursive", "statespace", "nonrobust"), ...) {
  newdata <- na.fail(newdata)
  if (n.ahead < 1L) stop("'n.ahead' must be at least 1")
  method <- match.arg(method)
  p <- object$order
  P <- seq_len(object$order)
  sd.pred <- sqrt(object$var.pred)
  if (is.null(object$x.intercept)) {
    xint <- 0
  } else {
    xint <- object$x.intercept
  } 
  if (p > 0) {
    # we do not need newdata if the model order is zero
    if (missing(newdata)) newdata <- object$x # use the observations to which the model was originally fitted if no other data are provided  
    st <- tsp(as.ts(newdata))[2L]
    dt <- deltat(newdata)
    xfreq <- frequency(newdata)
    tsp(newdata) <- NULL
    class(newdata) <- NULL
    newdata <- filtered.arrob(object, x = newdata, method = method, ...) # filtering
    n <- length(newdata)
    x <- c(newdata - object$x.mean, rep.int(0, n.ahead))
    for (i in seq_len(n.ahead)) {
      x[n + i] <- sum(object$ar * x[n + i - P]) + xint
    }
    pred <- x[n + seq_len(n.ahead)] + object$x.mean   
    if (se.fit) {
      psi <- if(n.ahead > 1) ARMAtoMA(ar = object$ar, lag.max = n.ahead - 1L) else NULL
      vars <- cumsum(c(1, psi^2))
      se <- (sd.pred * sqrt(vars))[seq_len(n.ahead)]
    }
  } else {
    pred <- rep.int(xint, n.ahead) + object$x.mean
    if (se.fit) se <- rep.int(sd.pred, n.ahead)
  }
  pred <- ts(pred, start = st + dt, frequency = xfreq)
  if (se.fit) {
    res <- list(pred = pred, se = ts(se, start = st + dt, frequency = xfreq))
  } else {
    res <- pred
  }
  return(res)
}

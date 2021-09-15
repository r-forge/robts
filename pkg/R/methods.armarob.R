print.armarob <- function (x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\n")
    cat("\n")
    cat("Estimated order of fitted ARMA model")
    cat("\n")
    cat(paste("AR part: ", x$order[1],sep=""))
    cat("\n")
    cat(paste("MA part: ", x$order[2],sep=""))
    cat("\n")
    cat("\n")
    cat("\n")
    cat("Residuals:")
    cat("\n")
    resid <- x$resid
    print(summary(resid), digits = digits)
    cat("\n")
    cat("\n")
    if (x$order[1] == 0) {
        cat("\n")
        cat("No AR Coefficients")
        cat("\n")
    }
    else {
        coefs <- matrix(x$ar,ncol=1)
        colnames(coefs) <- "AR coefficients"
        rownames(coefs) <- 1:length(length(x$ar))
        printCoefmat(coefs, digits = digits)
    }
    cat("\n")
    cat("\n")
    if (x$order[2] == 0) {
        cat("\n")
        cat("No MA Coefficients")
        cat("\n")
    }
    else {
        coefs <- matrix(x$ma,ncol=1)
        colnames(coefs) <- "MA coefficients"
        rownames(coefs) <- 1:length(length(x$ma))
        printCoefmat(coefs, digits = digits)
    }
    cat("\n")
    cat("\n")
    cat("Estimated location of time series:", format(signif(x$x.mean, 
        digits)))
    cat("\n")
    cat("Estimated variance of residuals:", format(signif(x$var.pred, 
        digits)))
    cat("\n")
    cat("\n")
    cat("\n")
    invisible(x)
}

predict.armarob <- function(object, newdata = object$x, n.ahead = 1, se.fit = TRUE, ...) {
  newdata <- na.fail(newdata)
  if (n.ahead < 1L) stop("'n.ahead' must be at least 1")
  sd.pred <- sqrt(object$var.pred)
  if (is.null(object$x.intercept)) {
    xint <- 0
  } else {
    xint <- object$x.intercept
  } 
  p <- object$order[1]
  q <- object$order[2]
  k <- max(p,q)
  if (k > 0) {
    # we do not need newdata if the model order is zero
    if (missing(newdata)) newdata <- object$x # use the observations to which the model was originally fitted if no other data are provided  
    st <- tsp(as.ts(newdata))[2L]
    dt <- deltat(newdata)
    xfreq <- frequency(newdata)
    tsp(newdata) <- NULL
    class(newdata) <- NULL
    if (p==0) {arv=0} else {arv <- object$ar}  
    if (q==0) {mavv=0} else {mav <- object$ma}
    filterednewdata <- filterrob.statespaceARMA(newdata,arv,mav,sd.pred^2) # filtering
    n <- length(newdata)
    x <- c(filterednewdata$filtered - object$x.mean, rep.int(0, n.ahead))
    pv <- max(1,p)
    qv <- max(1,q)
    u <- filterednewdata$residuals
    u <- c(u,rep(0,n.ahead))
    for (i in seq_len(n.ahead)) {
      x[n + i] <- as.numeric(t(x[n+i-(1:pv)])%*%arv) + as.numeric(t(u[n+i-(1:qv)])%*%mav)+ xint
    }
    pred <- x[n + seq_len(n.ahead)] + object$x.mean   
    if (se.fit) {
      psi <- if(n.ahead > 1) ARMAtoMA(ar = arv,ma = mav, lag.max = n.ahead - 1L) else NULL
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


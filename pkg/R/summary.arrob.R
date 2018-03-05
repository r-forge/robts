summary.arrob <- function(object,correlation=FALSE,symbolic.cor = FALSE,...) {
ans <- list()
class(ans) <- "summary.arrob"
if (object$order==0) ans$residuals <- object$resid else ans$residuals <- object$resid[(object$order+1):object$n.used]
ans$method <- object$method
ans$order <- object$order
ans$estsig <- object$var.pred
ans$mean <- object$x.mean
ans$coefvar <- object$asy.var.coef
ans$n <- object$n.used
ans$symbolic.cor <- symbolic.cor
if (!is.null(ans$coefvar)) {
    ster <- diag(ans$coefvar)
    teststat <- object$ar/sqrt(ster)
    ans$coefficients <- cbind(Estimate = object$ar, `Std. Error` = ster,`test stat.` = teststat, `Pr(>|t|)` = 2 * pnorm(abs(teststat), lower.tail = FALSE))
} else {ans$coefficients <- cbind(Estimate = object$ar, `Std. Error` = NA,`test stat.` = NA, `Pr(>|t|)` = NA)}
if (ans$order>0) {
    rowna <- character()
    for (i in 1:ans$order) rowna <- c(rowna,paste("AR",i))
    rownames(ans$coefficients) <- rowna
    }
if (correlation&(ans$order>0)) {
    ans$cor <- matrix(ncol=ans$order,nrow=ans$order)
    colnames(ans$cor) <- rowna
    rownames(ans$cor) <- rowna
    if (!is.null(object$asy.var.coef)) ans$cor[1:ans$order,1:ans$order] <- cov2cor(object$asy.var.coef)
} else ans$cor <- NULL
ans$corp <- correlation
return(ans)
}


print.summary.arrob <- function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) {
    cat("\n")
    cat("\n")
    cat(paste("Estimated order of fitted autoregressive model:",x$order))
    cat("\n")
    cat("\n")
    cat("\n")
    cat("Residuals:")
    cat("\n")
    resid <- x$residuals
    print(summary(resid), digits = digits, ...)
    cat("\n")
    cat("\n")
    if (x$order == 0) {
        cat("\n")
        cat("No Coefficients")
        cat("\n")
    }
    else {
        coefs <- x$coefficients
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\n")
    cat("\n")
    cat("Estimated location of time series:", format(signif(x$mean, 
        digits)))
    cat("\n")
    cat("Estimated scale of residuals:", format(signif(x$estsig, 
        digits)))
    cat("\n")
    cat("\n")
    correl <- x$cor
    if (x$corp){
    if (!is.null(correl)) {
        cat("Estimated correlation of parameters:")
        cat("\n")
        p <- NCOL(correl)
        if (p > 1L) {
            cat("\n")
            cat("Correlation of Coefficients:")
            cat("\n")
            if (is.logical(x$symbolic.cor) && x$symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }}
    cat("\n")
    invisible(x)
}



plot.acfrob <- function(x, ci = 0.95, ci.col = "blue", ci.type = "white", ...){
  stats:::plot.acf(x, ci = 0, ...) # plot a usual acf without confidence intervals
  with.ci <- ci > 0 && x$type != "covariance"
  if (with.ci){
    ci.type <- match.arg(ci.type)
    clim0 <- qnorm((1 + ci)/2) / sqrt(x$n.used) * x$are
    abline(h = c(-1, +1)*clim0, col = ci.col, lty = "dashed")
  }
  invisible()
}

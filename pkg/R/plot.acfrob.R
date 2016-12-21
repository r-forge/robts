plot.acfrob <- function(x, ci = 0.95, ci.col = "blue", ci.type = "white", ylim = NULL, ...){
  with.ci <- ci > 0 && x$type != "covariance"  
  if (with.ci){
    ci.type <- match.arg(ci.type)
    clim0 <- qnorm((1 + ci)/2) / sqrt(x$n.used) * x$are
  }
  if(is.null(ylim)) {
    ylim <- range(x$acf, na.rm = TRUE)
    if(with.ci && !is.na(clim0)) ylim <- range(-clim0, clim0, ylim)                  
  }
  plot(structure(x, class = "acf"), ci = 0, ylim = ylim, ...) # plot a usual acf without confidence intervals
  if (with.ci) abline(h = c(-1, +1)*clim0, col = ci.col, lty = "dashed")
  invisible()
}

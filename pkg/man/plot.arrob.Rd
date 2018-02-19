\encoding{utf8}
\name{plot.arrob}
\alias{plot.arrob}
\title{
Diagnostic Plots for a Fitted Autoregressive Model
}
\description{
Diagnostic plots for objects of class \code{"arrob"} as it is returned by a call to the function \code{\link{arrob}}.
}
\usage{
\method{plot}{arrob}(x, ask = TRUE, ci = 0.95, ...)
}
\arguments{
  \item{x}{
an object of class \code{"arrob"}. Usually the result of a call to \code{\link{arrob}}.
}
  \item{ask}{
logical value. If \code{TRUE} (and the R session is interactive) the user is asked for input, before a new figure is drawn (see \code{\link{devAskNewPage}}).
}
  \item{ci}{coverage probability for confidence intervals in some of the plots. Plotting of the confidence interval is suppressed if \code{ci} is zero or negative.}
  \item{...}{
further arguments are currently ignored. Only for compatibility with generic function.
}
}
\details{
Produces the following plots:
\itemize{
\item{Information criterion for different model orders. For robust fits obtained by \code{\link{arrob}} the plot shows a robust information criterion; for classical fits obtained by \code{\link{ar}} it shows differences of Akaikes Information criterion (AIC) with the AIC of the best model. The best model is marked by a vertical line. This plot is helpful for selection of the optimal order of an AR model.}
\item{Robust autocorrelation function of the residuals. The autocorrelation function is robustly estimated using the default settings of \code{\link{acfrob}}. If the confidence level \code{ci} is greater than zero, then confidence limits under the assumption of uncorrelated residuals are added to the plot, see \code{\link{plot.acfrob}}. This plot is helpful to check whether the residuals are uncorrelated.}
\item{Residuals plotted against time. If the confidence level \code{ci} is greater than zero, then confidence limits under the assumption of uncorrelated residuals with zero mean and constant variance are added to the plot. This plot is helful to check whether the residuals have mean zero and constant variance and to identify possible outliers.}
\item{Normal quantile-quantile plot of the residuals, see functions \code{\link{qqnorm}} and \code{\link{qqline}}. This plot is helpful to check whether the residuals are from a normal distribution.}
}

}
\author{
Alexander \enc{Dürre}{Duerre} and Tobias Liboschik
}

\seealso{
\code{\link{arrob}} for robust fitting of AR models, \code{\link{ar}} for classical fitting.

Other S3 methods \code{\link{residuals.arrob}}, \code{\link{fitted.arrob}}, \code{\link{filtered.arrob}} and \code{\link{predict.arrob}} for objects of class \code{"arrob"}. 
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 0.5), n = 100)
fit <- arrob(tss, method = "yw")
plot(fit)
}

\keyword{Autoregressive models}

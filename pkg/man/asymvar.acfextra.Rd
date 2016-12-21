\encoding{utf8}
\name{asymvar.acfextra}
\alias{asymvar.acfextra}

\title{
Calculation of the Long Run Variance Based on an AR Fit
}
\description{
Computes the long run variance, which is required for change point testing, by first fitting an AR-model and then extrapolating its theoretic autocovariance function.
}

\usage{
asymvar.acfextra(x, obs = c("untransformed", "ranks"), order.max = 2,
                 aic = FALSE)
}

\arguments{
	\item{x}{numeric vector or univariate time series object.}
	\item{obs}{character string indicating whether the long run variance of a cusum statistic (\code{"untransformed"}) or of a Wilcoxon statistic (\code{"ranks"}) is calculated.}
	\item{aic}{logical indicating whether the AR order should be estimated by robust AIC criteria considering orders up to \code{order.max} (if \code{aic = TRUE}) or should be equal to \code{order.max} (if \code{aic = FALSE}).}
	\item{order.max}{integer value determining the (maximal) order of the AR fit.}
}

\details{
Cusum-type change point tests require an estimation of the long run variance. One possibility is to fit an AR model and calculate its corresponding autocovariance function. The long run variance is then estimated by the Bartlett estimator with bandwidth \eqn{n-1}.
}

\value{
List containing the following named elements:
\item{lrv}{estimated long run variance}
\item{order}{estimated order of the AR model}
}

\seealso{
The long run variance can be also estimated by \code{\link{asymvar.window}} and \code{\link{asymvar.acf}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.acfextra(tss)
}

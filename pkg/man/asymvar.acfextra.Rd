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
asymvar.acfextra(x,obs=c("untransformed","ranks"),order.max=2,aic=FALSE,...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{obs}{character string indicating whether the long run variance of a cusum statistic or a Wilcoxon statistic is calculated.}
	\item{aic}{logical whether the order of the AR model is estimated by the AIC criterion or is considered as fixed.}
	\item{order.max}{If \code{aic=TRUE} this is the maximal considered AR model, otherwise it is the used order of the AR model.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
List containing the following values:
\item{lrv}{estimated long run variance}
\item{order}{estimated order of the AR model}
}

\details{
Cusum-type change point tests require an estimation of the long run variance. One possibility is to fit an AR model and calculate its corresponding autocovariance function. The long run variance is then estimated by the Bartlett estimator with bandwidth \eqn{n-1}.
}


\seealso{

The long run variance can be also estimated by \code{\link{asymvar.window}} and \code{\link{asymvar.acf}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.acfextra(tss)
}

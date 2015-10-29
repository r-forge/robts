\encoding{utf8}
\name{asymvar.acfextra}
\alias{asymvar.acfextra}

\title{
Calculation of the long run variance based on a an AR-fit
}
\description{
Computes the long run variance, which is required for change point testing, by first fitting an AR-model and then extrapoliting its theoretic autocovariancefunction.
}

\usage{
asymvar.acfextra(x,shiftcorrect=TRUE,obs=c("untransformed","ranks"),borderN=10,order.max=2,aic=FALSE,...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{shiftcorrect}{logical, if TRUE the estimation of the long run variance considers a change in location by first estimating the changepoint and then subtracting the estimated jump height from the second part of the time series}
	\item{obs}{character string indicating whether the lon run variance of a cusum statistic or a Wilcoxon statistic is calculated.}
	\item{borderN}{numeric, if \code{shiftcorrect == TRUE} the changepoint for the estimation of the log run variance is not searched below the first and last borderN observations.}
	\item{aic}{logical whether the order of the AR model is estimated by the aic criterium or is considered as fixed.}
	\item{order.max}{If \code{aic=TRUE} this is the maximal considered AR modell, otherwise it is the used order of the AR modell.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
List containing the folowing values:
\item{lrv}{estimated long run variance}
\item{order}{estimated order of the AR modell}
}

\details{
Cusum-type changepoint tests require an estimation of the long run variance. One possibility is to fit an AR model and calculates its corresponding autocovariance function. The long run variance is then estimated by the bartlett estimator with bandwith \eqn{n-1}.
}


\seealso{

The long run variance can be also estimated by \code{\link{asymvar.window}} and \code{\link{asymvar.acf}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.acfextra(tss)
}

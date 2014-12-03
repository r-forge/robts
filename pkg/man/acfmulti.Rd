\encoding{utf8}
\name{acfmulti}
\alias{acfmulti}

\title{
Autocorrelation function from multivariate correlation
}
\description{
Calculates the autocorrelation function of a time series using multivariate correlation estimators, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfmulti(timeseries, maxlag, method = "weightedMCD", ...)
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
	\item{method}{Multivariate correlation estimator to be used: \code{"weightedMCD"}, \code{"rawMCD"}, \code{"Stahel-Donoho"}, \code{"S"}, \code{"M"}, \code{"reweight"}, or \code{"Tyler"}, see \enc{Dürre}{Duerre} et al. (2014) for details.}
	\item{...}{Further arguments for the correlation estimators.}
}

\value{
Vector of length \code{maxlag} containing the autocorrelation function.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.
}
\author{
Alexander \enc{Dürre}{Duerre}
}


\seealso{
Other functions for calculating the autocorrelation function are \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, and \code{\link{trimmedcor}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfmulti(timeseries = ts, maxlag = 20)
}

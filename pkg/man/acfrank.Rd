\encoding{utf8}
\name{acfrank}
\alias{acfrank}

\title{
Autocorrelation function by ranks
}
\description{
Calculates the autocorrelation function of a time series using rank estimators, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfrank(timeseries, maxlag, method = "gaussian")
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
	\item{method}{Correlation estimator to be used: \code{"gaussian"}, \code{"spearman"}, \code{"kendall"}, \code{"quadrant"}, or \code{"masarotto"}, see \enc{Dürre}{Duerre} et al. (2014) for details.}
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
Other functions for calculating the autocorrelation function are \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrobfil}}, and \code{\link{trimmedcor}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrank(timeseries = ts, maxlag = 20)
}

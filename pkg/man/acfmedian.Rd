\encoding{utf8}
\name{acfmedian}
\alias{acfmedian}

\title{
Autocorrelation function based on mediancorrelation
}
\description{
Calculates the autocorrelation function of a time series based on mediancorrelation, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfmedian(timeseries, maxlag)
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
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
Other functions for calculating the autocorrelation function are \code{\link{acfGK}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, and \code{\link{acftrim}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfmedian(timeseries = ts, maxlag = 20)
}

\encoding{utf8}
\name{robspec.acf}
\alias{robspec.acf}

\title{
Spectral density from autocorrelation function
}
\description{
Calculates the spectral density of a time series from its given autocorrelation function calculated robustly.
}
\usage{
robspec.acf(acf, tmax = length(acf) + 2)
}

\arguments{
	\item{acf}{Autocorrelation function of a time series calculated robustly.}
	\item{tmax}{Length of the time series.}
}

\value{
	An object of class \code{"spec"}, see \code{\link{spectrum}}, a named list with the essential components:
	\item{freq}{Vector of frequencies at which the spectral desity is evaluated: \code{floor(tmax / 2) + 1} points from 0 to 1 / 2.}
	\item{spec}{Numeric vector containing the robust spectral density of the time series evaluated at \code{freq}.}
}

\author{
Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acorf <- acf(tss, plot = FALSE)$acf[-1]
robspec.acf(acorf, tmax = 100)
}

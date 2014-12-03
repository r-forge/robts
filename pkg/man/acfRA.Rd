\encoding{utf8}
\name{acfRA}
\alias{acfRA}

\title{
Autocorrelation function by RA approach
}
\description{
Calculates the autocorrelation function of a time series using a RA approach, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfRA(timeseries, maxlag, Psi = "Huber", meanvalue = median, scattervalue = mad, ...)
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
	\item{Psi}{Psi function to be used: \code{"Huber"} or \code{"Tukey"}, see \enc{Dürre}{Duerre} et al. (2014) for details.}
	\item{meanvalue}{Mean estimator to be used.}
	\item{scattervalue}{Scale estimator to be used.}
	\item{...}{Robustness parameter for \code{Psi} functions.}
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
Other functions for calculating the autocorrelation function are \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, and \code{\link{trimmedcor}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfRA(timeseries = ts, maxlag = 20)
}

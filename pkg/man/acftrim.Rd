\encoding{utf8}
\name{acftrim}
\alias{acftrim}

\title{
Autocorrelation function by trimmimg
}
\description{
Calculates the autocorrelation function of a time series by trimming it, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acftrim(timeseries, maxlag, alpha = 0.1)
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
	\item{alpha}{Part of largest and smallest observations to be trimmed.}
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
Other functions for calculating the autocorrelation function are \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, and \code{\link{acfrobfil}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acftrim(timeseries = ts, maxlag = 20)
}

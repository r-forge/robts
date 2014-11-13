\encoding{utf8}
\name{acfGK}
\alias{acfGK}

\title{
Autocorrelation function by GK method
}
\description{
Calculates the autocorrelation function of a time series using the GK method: estimation based on variances, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfGK(timeseries, maxlag, method = c("Qn", "MAD", "Tau"))
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
	\item{method}{Variance estimator to be used, see \enc{Dürre}{Duerre} et al. (2014) for details.}
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
Other functions for calculating the autocorrelation function are \code{\link{acfmedian}}, \code{\link{acfyohai}}, \code{\link{Autokor}}, \code{\link{multi}}, \code{\link{RAacf}}, \code{\link{rankacf}}, and \code{\link{trimmedcor}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfGK(timeseries = ts, maxlag = 40)
}

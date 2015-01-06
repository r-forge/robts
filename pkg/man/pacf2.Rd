\encoding{utf8}
\name{pacf2}
\alias{pacf2}

\title{Partial autocorrelations}

\description{
Calculates the partial autocorrelations of a time series using rank estimators, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
pacf2(timeseries, maxlag, method = "spearman")
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{maxlag}{Maximum lag of interest.}
	\item{method}{Correlation estimator to be used: \code{"spearman"}, \code{"gaussian"}, \code{"kendall"}, \code{"quadrant"}, or \code{"masarotto"}, see \enc{Dürre}{Duerre} et al. (2014) for details.}
}

\value{
Vector of length \code{maxlag} containing the partial autocorrelations.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.
}
\author{
Alexander \enc{Dürre}{Duerre}
}


\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
pacf2(timeseries = ts, maxlag = 20)
}

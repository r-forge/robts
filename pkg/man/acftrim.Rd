\encoding{utf8}
\name{acftrim}
\alias{acftrim}

\title{
Robust Autocorrelation Estimation Based on Trimming
}
\description{
Robustly estimates the autocorrelation function of a time series based on the trimmed timeseries. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acftrim(timeseries,maxlag,alpha=0.1)
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
	\item{alpha}{numeric indicating the amount of trimmed values, see details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

The function basically first trimms the smallest and largest observations of the time series and then calculates the usual sample acf. This procedure was proposed by Chan and Wei (1992).

\code{alpha} determines the quantiles (\eqn{\alpha} and \eqn{1-\alpha}) at which the data is trimmed. Therefore \eqn{2\alpha} of the data is trimmed, see \enc{Dürre}{Duerre} et al. (2014) for details.

Since the estimator does not measure the usual autocorrelation, a nonlinear transformation based on simulations is applied, which guarantees consistency under gaussian processes.


}

\references{

Chan, W., and Wei, W. (1992): A comparison of some estimators of time series autocorrelations, \emph{Computational statistics & data analysis}, vol 14, .149--163.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrobfil}}, \code{\link{acfrank}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acftrim(tss,7)
}

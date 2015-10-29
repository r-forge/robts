\encoding{utf8}
\name{acfmedian}
\alias{acfmedian}

\title{
Robust autocorrelation estimation based on median correlation
}
\description{
Robustly estimates the autocorrelation function of a time series based on median correlation. See Chakhchoukh (2010) for details.
}

\usage{
acfmedian(timeseries,maxlag)
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

This function estimates the autocorrelation function based on the median correlation, which was proposed by Chakhchoukh (2010). The acf is estimated for each lag individually. Let X denote the original timeseries and Y the lagged one, then the autocorrelation is estimated by
\deqn{median(XY)/median(X^2).}
This value does not measure the usual autocorrelation. Therefore a nonlinear transformation based on a Monte Carlo simulation is applied, which is only valid under gaussian processes.
}  

\references{
Chakhchoukh, Y. (2010): A new robust estimation method for ARMA models, \emph{IEEE Transactions on Signal Processing}, vol. 58, 3512--3522.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfmedian(tss,7)
}

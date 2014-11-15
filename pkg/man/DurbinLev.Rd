\encoding{utf8}
\name{DurbinLev}
\alias{DurbinLev}

\title{
Partial autocorrelation function by Durbin-Levinson algorithm
}
\description{
Calculates the partial autocorrelation function using the Durbin-Levinson, see Maronna et al. (2006), chapter 8.2.1 for details.
}
\usage{
DurbinLev(ts, acf, robvar = TRUE)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{acf}{Autocorrelation function calculated from \code{ts}.}
	\item{robvar}{Logical, whether the variance of \code{ts} is estimated robustly. If \code{TRUE}, \code{mad} is used, else \code{sd}.}
}

\value{
Named list of 2:
\item{phis}{List of length of \code{acf} containing the coefficient vectors.}
\item{nus}{Numeric vector of same length + 1 containing the mean squared errors.}
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Tobias Liboschik and Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acor <- acfGK(timeseries = tss, maxlag = 6)
DurbinLev(ts = tss, acf = acor)
}

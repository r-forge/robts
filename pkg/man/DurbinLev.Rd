\encoding{utf8}
\name{DurbinLev}
\alias{DurbinLev}

\title{
Partial autocorrelation function by Durbin-Levinson algorithm
}
\description{
Calculates the partial autocorrelation function using the Durbin-Levinson algorithm, see Maronna et al. (2006), chapter 8.2.1 for details.
}
\usage{
DurbinLev(acf)
}

\arguments{
	\item{acf}{Autocorrelation function.}
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
DurbinLev(acf = acor)
}
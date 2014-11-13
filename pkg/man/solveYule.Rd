\encoding{utf8}
\name{solveYule}
\alias{solveYule}

\title{
AR model coefficients by Yule-Walker equations
}
\description{
Calculates the regression coefficients of an AR model by solving the Yule-Walker equations, see Maronna et al. (2006), chapter 8.2 for details.
}
\usage{
solveYule(ts, acf, p, robvar = TRUE)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{acf}{Autocorrelation function calculated from \code{ts}.}
	\item{p}{Order of the AR model.}
	\item{robvar}{Logical, whether the variance of \code{ts} is estimated robustly. If \code{TRUE}, \code{mad} is used, else \code{sd}.}
}

\value{
Vector of length \code{p} containing the regression coefficients of the AR(p) model.
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acor <- acfGK(timeseries = tss, maxlag = 6)
solveYule(ts = tss, acf = acor, p = 6)
}

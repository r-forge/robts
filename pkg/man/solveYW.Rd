\encoding{utf8}
\name{solveYW}
\alias{solveYW}

\title{
AR model coefficients by Yule-Walker equations
}
\description{
Calculates the regression coefficients of an AR model by solving the Yule-Walker equations, see Maronna et al. (2006), chapter 8.2 for details.
}
\usage{
solveYW(acf, p)
}

\arguments{
	\item{acf}{Autocorrelation function calculated from \code{ts}.}
	\item{p}{Order of the AR model.}
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
acor <- as.numeric(acfrob(tss, lag.max = 6, plot = FALSE)$acf)
solveYW(acf = acor, p = 6)
}
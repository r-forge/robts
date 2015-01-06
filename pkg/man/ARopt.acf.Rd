\encoding{utf8}
\name{ARopt.acf}
\alias{ARopt.acf}

\title{
Optimal order of AR model calculated via the autocorrelation function
}
\description{
Calculates the regression coefficients of an AR(p) model optimized with respect to p via the autocorrelation function using the Durbin-Levinson algorithm.
}

\details{
The regression coefficients and a robust version of Akaikes information criterion are calculated for every p. The solution with the smallest AIC is returned.

See also Maronna et al. (2006), chapter 8.2 and 8.6.6.
}

\usage{
ARopt.acf(tss, pmax = NULL, acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"))
}

\arguments{
	\item{tss}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{pmax}{Maximum order of AR model to consider. If \code{NULL}, the default minimum of \eqn{(T - 1) / 4} and \eqn{10 * lg(T)} is used.}
	\item{acf.fun}{Character value naming the function to calculate the autocorrelation function. See \code{"\link{acfrob}"} for details. Do not use \code{"acfrobfil"}.}
}

\value{
Numeric vector containing the regression coefficients of the optimal AR(p) model.
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Jonathan Rathjens and Tobias Liboschik
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
ARopt.acf(tss)
}

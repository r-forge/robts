\encoding{utf8}
\name{robspec}
\alias{robspec}

\title{
Robust spectrum
}
\description{
Calculates an estimate of the spectral density of a time series robustly. See Maronna et al. (2006), chapter 8.14.
}

\details{
The AR(p) model coefficients and the respective optimal p ist computed by \code{\link{ARopt.acf}}. The smoothed periodogram is calculated from the filtered residuals of this model.
}

\usage{
robspec(tss, psifunc = smoothpsi, acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acftrim"), spans = 8)
}

\arguments{
	\item{tss}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{psifunc}{Psi function for filtering of residuals, currently no alternative.}
	\item{acf.fun}{Character value naming the function to calculate the autocorrelation function. See \code{"\link{acfrob}"} for details. Do not use \code{"acfrobfil"}.}
	\item{spans}{Window width for smoothing kernel, has to be an even number.}	
}

\value{
	An object of class \code{"spec"}, see \code{\link{spectrum}}, a named list with the essential components:
	\item{freq}{Vector of frequencies at which the spectral desity is evaluated: \code{floor(tmax / 2) + 1} points from 0 to 1 / 2 with \code{tmax = length(tss) - p}.}
	\item{spec}{Numeric vector containing the robust spectral density of the time series evaluated at \code{freq}.}
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 50)
robspec(tss)
}

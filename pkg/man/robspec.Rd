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
The AR(p) model coefficients and the respective optimal p ist computed by \code{\link{lmrobARopt}}. The smoothed periodogram is calculated from the residuals of this model filtered by \code{\link{ts.robfilter}}. See the respective manual.

The resulting spectrum is evaluated at \code{length(ts)} points from - 1 / 2 to + 1 / 2.
}

\note{
The smoothing function returns a vector of the length put in. Consider trimming the ends of the final result.
}

\usage{
robspec(ts,
	smooth.type = "runmean", M = 1, taper = "cosine",
	method = "MM", singular.ok = FALSE, init = NULL, 
	psifunc = smoothpsi)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s. Its length should be not less than 7.}
	\item{smooth.type}{Smoothing function: default \code{"runmean"} for running mean, \code{"none"} for no smoothing.}
	\item{M}{Half window width for smoothing.}
	\item{taper}{Data taper: default \code{"cosine"}, or \code{"none"}.}
	\item{method}{String specifying the estimator-chain, to be passed to \code{\link{lmrob}}, see there for details.}
	\item{singular.ok}{Logical: if \code{FALSE}, a singular fit is an error, to be passed to \code{\link{lmrob}}, see there for details.}
	\item{init}{Initial estimate to be passed to \code{\link{lmrob}}, see there for details.}
	\item{psifunc}{Psi function to be passed to \code{\link{ts.robfilter}}, currently no alternative.}
}

\value{
Vector of length \code{tmax = length(ts)} (+ 1 if \code{tmax} is even) containing the estimated spectral density evaluated at \code{tmax} (+ 1) points from - 1 / 2 to + 1 / 2.
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

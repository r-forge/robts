\encoding{utf8}
\name{acfGK}
\alias{acfGK}

\title{
Robust Autocorrelation Function Estimation by GK Method
}
\description{
Robustly estimates the autocorrelation function of a time series using the Gnanadesikan-Kettenring (GK) method: estimation based on variances, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfGK(timeseries, maxlag, method = "Qn", ...)
}

\arguments{
	\item{timeseries}{a univariate numeric vector or time series object. Must not contain missing values.}
	\item{maxlag}{maximum lag at which to calculate the acf.}
	\item{method}{character string giving the type of variance estimator to be used. Possible values are \code{"Qn"}, \code{"MAD"} or \code{"Tau"}. See \enc{Dürre}{Duerre} et al. (2014) for details.}
	\item{...}{additional arguments passed to the function \code{\link[robustbase]{scaleTau2}} in package \code{robustbase}.}
}

\value{
Vector of length \code{maxlag} containing the robustly estimated autocorrelation function.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.
}
\author{
Alexander \enc{Dürre}{Duerre}
}


\seealso{
Other functions for robustly estimating the autocorrelation function are \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}} and \code{\link{acftrim}}. Classical, nonrobust estimation is provided by the function \code{\link{acf}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfGK(timeseries = ts, maxlag = 20)
}

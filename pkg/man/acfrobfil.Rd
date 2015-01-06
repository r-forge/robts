\encoding{utf8}
\name{acfrobfil}
\alias{acfrobfil}

\title{
Autocorrelation function from robustly filtered values
}
\description{
Calculates the autocorrelation function of a time series using robustly filtered values, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfrobfil(timeseries, p, maxlag, psifunc = smoothpsi)
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{p}{Order of AR fit to be used.}
	\item{maxlag}{Maximum lag of interest.}
	\item{psifunc}{Psi function to be used.}
}

\value{
	List of 8:
	\item{1}{Autocorrelation function using the empirical acf of the robustly filtered timeseries.}
	\item{2}{Autocorrelation function using the estimated partial autocorrelations.}
	\item{3}{Autocorrelation function using the estimated partial autocorrelations after Mottonen et al. (1999).}
	\item{4}{AIC value with respect to 1.}
	\item{5}{AIC value with respect to 2.}
	\item{6}{AIC value with respect to 3.}
	\item{7}{Order of AR fit used for robust filtering.}
	\item{8}{Robustly filtered time series.}
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.

Mottonen, J., Koivunen, V., and Oja, H. (1999): Robust autocovariance estimation based on sign and rank correlation coefficients, In: Proceedings of the IEEE Signal Processing Workshop on Higher-Order Statistics, 187–190.
}
\author{
Alexander \enc{Dürre}{Duerre}
}


\seealso{
Other functions for calculating the autocorrelation function are \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, and \code{\link{acftrim}}.
}

\examples{
set.seed(1066)
ts <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrobfil(timeseries = ts, maxlag = 20)
}

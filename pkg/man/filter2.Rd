\encoding{utf8}
\name{filter2}
\alias{filter2}

\title{
Robustly filtered time series using AR model coefficients
}
\description{
Filters a given time series robustly using the coefficients of a fitted AR model, see Maronna et al. (2006), chapter 8.6.2 for details.
}
\usage{
filter2(timeseries, phi, psifunc = smoothpsi)
}

\arguments{
	\item{timeseries}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{phi}{Coefficients of the AR model fitted to \code{timeseries}.}
	\item{psifunc}{Psi function, currently no alternative.}
}

\value{
	Named list of 2:
	\item{filtered.ts}{The robustly filtered version of \code{timeseries}.}
	\item{residuals}{The residuals of the fitted AR model.}
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Alexander \enc{Dürre}{Duerre} and Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acorf <- as.numeric(acfrob(tss, lag.max = 3, plot = FALSE)$acf)
phis <- solveYule(acorf, p = 3)
filter2(timeseries = tss, phi = phis)
}

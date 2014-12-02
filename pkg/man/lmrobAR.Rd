\encoding{utf8}
\name{lmrobAR}
\alias{lmrobAR}

\title{
AR model coefficients from robustly fitted linear model
}
\description{
Calculates the regression coefficients of an AR model by robustly fitting a linear model, see Maronna et al. (2006), chapter 8.2 for details. The model fit is performed by \code{\link{lmrob}}, see Rousseeuw et al. (2014).
}
\usage{
lmrobAR(ts, p, interc = TRUE, singular.ok = FALSE, ...)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{p}{Order of the AR model.}
	\item{interc}{Logical, whether the model contains an intercept, i. e., the time series is not centered to a zero mean.}
	\item{singular.ok}{Logical: if \code{FALSE}, a singular fit is an error, to be passed to \code{\link{lmrob}}, see there for details.}
	\item{...}{Further arguments to be passed to \code{\link{lmrob}}, see there for details.}
}

\value{
Named list of 2:
\item{coefficients}{Vector of length \code{p} (+ 1 if \code{interc = TRUE}) containing the regression coefficients of the AR(p) model (with the intercept in the last place).}
\item{model}{Object of class \code{lmrob}: the fitted model, see \code{\link{lmrob}} for details.}
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.

Rousseeuw, P., Croux, C., Todorov, V., Ruckstuhl, A., Salibian-Barrera, M., Verbeke, T., Koller, M., and Maechler, M. (2014): robustbase: Basic Robust Statistics. R package version 0.91-1. URL \url{http://CRAN.R-project.org/package=robustbase}.
}
\author{
Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
lmrobAR(tss, 6)
}

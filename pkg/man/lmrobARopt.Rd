\encoding{utf8}
\name{lmrobARopt}
\alias{lmrobARopt}

\title{
AR model coefficients from robustly fitted AIC-optimal linear model
}
\description{
Calculates the regression coefficients of an AR(p) model by robustly fitting a linear model optimized with respect to p by the AIC.
}

\details{
The function \code{\link{lmrobAR}} is applied for every p, with a robust version of Akaikes information criterion calculated each, at last returning the solution with the smallest AIC. Models without convergence are not considered.

See also Maronna et al. (2006), chapter 8.2 and 8.6.6.

The model fit is performed by \code{\link{lmrob}}, see Rousseeuw et al. (2014).
}

\usage{
lmrobARopt(ts, interc = TRUE, singular.ok = FALSE, pmax = NULL, ...)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s. Its length should be not less than 7.}
	\item{interc}{Logical, whether the model contains an intercept, i. e., the time series is not centered to a zero mean.}
	\item{singular.ok}{Logical: if \code{FALSE}, a singular fit is an error, to be passed to \code{\link{lmrob}}, see there for details.}
	\item{pmax}{Maximum order of AR model to consider. If \code{NULL}, the default minimum of \eqn{(T - 2) / 4} and \eqn{10 * lg(T)} is used.}
	\item{...}{Further arguments to be passed to \code{\link{lmrob}}, see there for details.}
}

\value{
Named list of 3:
\item{coefficients}{Vector of length p (+ 1 if \code{interc = TRUE}) containing the regression coefficients of the optimal AR(p) model (with the intercept in the last place).}
\item{model}{Object of class \code{lmrob}: the optimal fitted model, see \code{\link{lmrob}} for details.}
\item{p.optimal}{The optimal p.}
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
lmrobARopt(tss)
}

\encoding{utf8}
\name{arrob}
\alias{arrob}
\alias{arrob.yw}
\alias{arrob.regression}
\alias{arrob.filter}

\title{
Robust Fitting of Autoregressive Models
}
\description{
	Robust fit of an autoregressive model to a time series.  Its order is optimized with respect to a robustified AIC or given.
}
\usage{
arrob(x, aic = TRUE, order.max = NULL,
      method = c("yw", "regression", "gm", "filter"),
      na.action = na.fail, series = deparse(substitute(x)), ...)

arrob.yw(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p,
      na.action = na.fail, series = deparse(substitute(x)),
      acf.approach = c("GK", "median", "multi", "partrank", "RA", "rank",
      "filter", "trim", "bireg"), locfn = median, scalefn = Qn, ...)
      
arrob.regression(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p,
      na.action = na.fail, series = deparse(substitute(x)),
      intercept = TRUE, scalefn = Qn, ...)
      
arrob.filter(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p,
      na.action = na.fail, series = deparse(substitute(x)),
      psi.l = 2, psi.0 = 3)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{aic}{logical indicating whether the AR order should be estimated by robust AIC criteria considering orders up to \code{order.max} (if \code{aic = TRUE}) or should be equal to \code{order.max} (if \code{aic = FALSE}).}
	\item{order.max}{integer value determining the (maximal) order of the AR fit.}
	\item{method}{character string naming the estimation method to be used, see Details.}
	\item{na.action}{function to be called to handle missing values.}
	\item{series}{the time series name.}
	\item{aicpenalty}{function of the model order defining the penalty term of the model selection criteria. The default results in a robust AIC.}
	\item{...}{further arguments to be passed to internal function \code{arrob.yw}, \code{arrob.regression}, \code{arrob.filter} or \code{\link{arrob.gm}}. See below or on the respective help page.}
	\item{acf.approach}{character string naming the function to calculate the autocorrelation function.  See \code{\link{acfrob}} for details.  It is not recommended to use \code{acf.approach = "filter"} since one can get an AR fit more directly by \code{method = "filter"}.}
	\item{locfn}{function which calculates a location estimator. Its argument must be a vector of the data.}
	\item{scalefn}{function which calculates a scale estimator. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
  \item{intercept}{logical whether the regression model includes an intercept.}
	\item{psi.l}{numeric value determining the used psi function, see the help page of \code{\link{ARfilter}} for details.}
	\item{psi.0}{numeric value determining the used psi function, see the help page of \code{\link{ARfilter}} for details.}
}

\details{
There are many possibilities for robust estimation of AR models. The argument \code{method} specifies which of the following methods is to be used, where \code{method = "yw"} is the default.

\describe{
\item{\code{"yw"}}{One approach is to use a roust estimation of the autocorrelation function. The corresponding AR coefficients can be then either computed by solving the Yule-Walker equations \code{method = "yule-walker"} employing the Durbin-Levinson algorithm. The acf estimator can be determined by \code{acf.fun}. See the help page of \code{\link{acfrob}} for an overview of implemented procedures for estimation of the acf.}

\item{\code{"regression"}}{The AR coefficients are computed by a robust regression. Here the dependent variable is the time series \code{x} and the independent variables are the up to \code{order.max} lagged time series. The regression is done by MM estimation with the function \code{\link{lmrob}}, see Rousseeuw et al. (2014).}

\item{\code{"filter"}}{The regression coefficients are computed via robust filtering as described in chapter 8.6 of Maronna et al. (2006). See the help page of \code{\link{ARfilter}} for details.}

\item{\code{"gm"}}{Generalized M estimates are used for obtaining an AR fit as described in Maronna et al. (2006). See the help page of \code{\link{arrob.gm}} for details and further arguments which are passed to that function by the \code{...} argument.}
}

If \code{aic = TRUE}, AR models for every order up to \code{order.max} are estimated and the one with the smallest model selection criterion returned. The criterion is the sum of the logarithm of a robust variance estimation of the residuals and a penalty term \eqn{pen(p)} divided by \eqn{n-p} depending on the AR order \eqn{p}, i.e.,
\deqn{log(\hat{\sigma}^2)+\frac{pen(p)}{n-p}.}{log(sigma_hat^2) + (pen(p))/(n-p)}
The penalty term can be set by the argument \code{aicpenalty}.
}

\value{
An object of class \code{"ar"}, see \code{\link{ar}} for details.

Note that list element \code{aic} gives the value of the robust information criterion and not its difference with the lowest information criterion of all considered models, as it is returned by the function \code{\link{ar}}.
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.

Rousseeuw, P., Croux, C., Todorov, V., Ruckstuhl, A., Salibian-Barrera, M., Verbeke, T., Koller, M. and Maechler, M. (2014): robustbase: Basic Robust Statistics. R package version 0.91-1. URL \url{http://cran.r-project.org/package=robustbase}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
Classical, nonrobust fitting is provided by the function \code{\link{ar}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
arrob(tss, method = "yw")
arrob(tss, method = "regression")
arrob(tss, method = "filter")
arrob(tss, method = "gm")
}

\keyword{Autoregressive models}

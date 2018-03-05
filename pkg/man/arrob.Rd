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
	Robust fit of an autoregressive model to a time series. Its order is either optimized with respect to a robustified AIC or specified by the user.
}
\usage{
arrob(x, aic = TRUE, order.max, method = c("yw", "regression", "gm", "filter"),
      na.action = na.fail, series = deparse(substitute(x)),asyvar=FALSE,bootpar=list(num=100,blockl=floor(2*length(x)^(1/2))), ...)

arrob.yw(x, order.max, aic = TRUE, aicpenalty = function(p) 2*p,
      na.action = na.fail, series = deparse(substitute(x)),
      acf.approach = c("GK", "median", "multi", "partrank", "RA", "rank",
      "filter", "trim", "bireg"), locfn = median, scalefn = Qn, ...)
      
arrob.regression(x, order.max, aic = TRUE, aicpenalty = function(p) 2*p,
      na.action = na.fail, series = deparse(substitute(x)),
      intercept = TRUE, scalefn = Qn, ...)
      
arrob.filter(x, order.max, aic = TRUE, aicpenalty = function(p) 2*p,
      na.action = na.fail, series = deparse(substitute(x)),
      psi.l = 2, psi.0 = 3)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{aic}{logical indicating whether the AR order should be estimated by robust AIC criteria considering orders up to \code{order.max} (if \code{aic = TRUE}) or should be equal to \code{order.max} (if \code{aic = FALSE}).}
	\item{order.max}{integer value determining the (maximal) order of the AR fit. If missing, this value is chosen to be \code{floor(min(c((n - 1) / 4, 10 * log(n, base = 10))))}, where \code{n} is the number of observations.}
	\item{method}{character string naming the estimation method to be used, see Details.}
	\item{na.action}{function to be called to handle missing values. Default is \code{\link{na.fail}}. Other reasonable possibilities are \code{\link{na.omit}} (omits all missing values at the beginning and the end of the time series and fails if there are more missing values), \code{\link{na.contiguous}} (uses longest contiguous stretch of non-missing observations) and \code{\link{na.extremify}} (experimental, sets missing values to very extreme values in the expectation that they have only little influence on the outcome of robust estimation methods).}
	\item{series}{the time series name.}
    \item{asyvar}{logical indicating whether the covariance matrix of the fitted AR coefficients should be estimated by a block bootstrap. The default is \code{FALSE} since it can take some time.}
    \item{bootpar}{list containing parameters of the blockbootstrap which is used to estimate the variance of the fitted AR coefficients, see Details.}
	\item{aicpenalty}{function of the model order defining the penalty term of the model selection criteria. The default results in a robust AIC.}
	\item{...}{further arguments to be passed to internal function \code{arrob.yw}, \code{arrob.regression}, \code{arrob.filter} or \code{\link{arrob.gm}}. See below or on the respective help page.}
	\item{acf.approach}{character string naming the function to calculate the autocorrelation function.  See \code{\link{acfrob}} for details.  It is not recommended to use \code{acf.approach = "filter"} since one can get an AR fit more directly by \code{method = "filter"}.}
	\item{locfn}{function which calculates a location estimator. Its argument must be a vector of the data.}
	\item{scalefn}{function which calculates a scale estimator. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
  \item{intercept}{logical whether the regression model includes an intercept.}
	\item{psi.l}{numeric value determining the psi function, see the help page of \code{\link{ARfilter}} for details.}
	\item{psi.0}{numeric value determining the psi function, see the help page of \code{\link{ARfilter}} for details.}
}

\details{
There are many possibilities for robust estimation of AR models. The argument \code{method} specifies which of the following methods is to be used, where \code{method = "yw"} is the default.

\describe{
\item{\code{"yw"}}{One approach is to use a roust estimation of the autocorrelation function. The corresponding AR coefficients can be then either computed by solving the Yule-Walker equations \code{method = "yule-walker"} employing the Durbin-Levinson algorithm. The acf estimator can be determined by \code{acf.fun}. See the help page of \code{\link{acfrob}} for an overview of implemented procedures for estimation of the acf.}

\item{\code{"regression"}}{The AR coefficients are computed by a robust regression. Here the dependent variable is the time series \code{x} and the independent variables are the up to \code{order.max} lagged time series. The regression is done by MM estimation with the function \code{\link{lmrob}}, see Rousseeuw et al. (2014).}

\item{\code{"filter"}}{The regression coefficients are computed via robust filtering as described in Chapter 8.6 of Maronna et al. (2006). See the help page of \code{\link{ARfilter}} for details.}

\item{\code{"gm"}}{Generalized M estimates are used for obtaining an AR fit as described in Maronna et al. (2006). See the help page of \code{\link{arrob.gm}} for details and further arguments which are passed to that function by the \code{...} argument.}
}

If \code{aic = TRUE}, AR models for every order up to \code{order.max} are estimated and the one with the smallest model selection criterion returned. The criterion is the sum of the logarithm of a robust variance estimation of the residuals and a penalty term \eqn{pen(p)} divided by \eqn{n-p} depending on the AR order \eqn{p}, i.e.,
\deqn{log(\hat{\sigma}^2)+\frac{pen(p)}{n-p}.}{log(sigma_hat^2) + (pen(p))/(n-p)}
The penalty term can be set by the argument \code{aicpenalty}.

If \code{asyvar=TRUE} the variance of the fitted AR coefficients is estimated by a moving block bootstrap, see Lahiri (1999) for a detailed description of this method. The blocklength can be set by the \code{blockl} argument and thr number of bootstrap repetitions by \code{num}. Note that there are no theoretical justifications that this estimation is consistent for robust methods.
}

\value{
Object of classes \code{"arrob"} and \code{"ar"}. This is a list which includes all elements of an object of class \code{"ar"} (see \code{\link{ar}} for details) plus the following additional element:

\item{x}{the original time series \code{x}.}

Note that list element \code{aic} gives the value of the robust information criterion and not its difference with the lowest information criterion of all considered models, as it is returned by the function \code{\link{ar}}.
}

\references{
Lahiri, S.N. (1999): Theoretical Comparisons of Block Bootstrap Methods, \emph{The annals of Statistics}, vol. 27, 386--404. 

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8, \doi{10.1002/0470010940}.

Rousseeuw, P., Croux, C., Todorov, V., Ruckstuhl, A., Salibian-Barrera, M., Verbeke, T., Koller, M. and Maechler, M. (2014): robustbase: Basic Robust Statistics. R package version 0.91-1. URL \url{http://cran.r-project.org/package=robustbase}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
Classical, nonrobust fitting is provided by the function \code{\link{ar}}.

S3 methods \code{\link{residuals.arrob}}, \code{\link{fitted.arrob}}, \code{\link{filtered.arrob}} and \code{\link{predict.arrob}} for objects of class \code{"arrob"}. 
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 0.5), n = 100)
arrob(tss, method = "yw")
arrob(tss, method = "regression")
arrob(tss, method = "filter")
arrob(tss, method = "gm")
}

\keyword{Autoregressive models}

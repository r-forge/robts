\encoding{utf8}
\name{arrob}
\alias{arrob}

\title{
Robust Estimation of Autoregressive Models
}
\description{
	Robust fit of an autogression model to a time series.  Its order is optimized with respect to a robustified AIC or given.
}
\usage{
arrob(x, aic = TRUE, order.max = NULL, method =
    c("yule-walker", "durbin-levinson", "robustregression", "filter", "gm"),
    na.action = na.fail, series = deparse(substitute(x)), ...,
    acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank",
        "acfRA", "acfrank", "acftrim"), 
    aicpenalty = function(p) {return(2 * p)})
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{aic}{logical whether the AR model order is to be optimized.}
	\item{order.max}{numeric value: if \code{aic == TRUE}, the maximum AR model order considered; if \code{aic == FALSE}, the given order.}
	\item{method}{character string naming the computation method used, see details.}
	\item{na.action}{function to be called to handle missing values.}
	\item{series}{the time series name.}
	\item{...}{further arguments to be passed to internal functions called.}
	\item{acf.fun}{character string naming the function to calculate the autocorrelation function.  See \code{\link{acfrob}} for details.  Do not use \code{"acfrobfil"} since one can get an AR-fit more directly by the method \code{filter}.}
	\item{aicpenalty}{function of the model order defining the penalty term of the model selection criteria. The default results in a robust AIC.}
}

\details{

There are many possibilities for robust estimation of AR models.

One approach is to use a roust estimation of the autocorrelation function. The corresponding AR coefficients can be then either computed by solving the Yule-Walker equations \code{method == "yule-walker"} or by using the Durbin-Levinson algorithm \code{method == "durbin-levinson"}. The acf estimator can be determined by \code{acf.fun}. See the help of \code{\link{acfrob}} for an overview of implemented procedures.

For \code{method == "robustregression"}, the AR coefficients are computed by a robust regression. Here the dependent variable is the timeseries and the \code{p} independent variables are the lagged timeserieses. The regression is done by MM estimator \code{\link{lmrob}}, see Rousseeuw et al. (2014).

For \code{method == "gm"}, generalized M estimates are used for the AR fit as described in Maronna et al. (2006), see the help page of \code{\link{bestAR}} for details and further arguments.

For \code{method == "filter"}, the regression coefficients are computed via robust filtering as described in the chapter 8.6 of Maronna et al. (2006), see the help page of \code{\link{ARfilter}} for details.

If \code{aic == "TRUE"} AR models for every order up to \code{order.max} are estimated and the one with the smallest model selection criteria returned. The criterion is the sum of the logarithm of a robust variance estimation of the residuals and a penalty term devidided by \code{n-p} depending on the AR order \code{p}
\deqn{crit=log(\hat{\sigma}^2)+\frac{pen(p)}{n-p}.}
The penalty term can be set by \code{aicpenalty}.
}
\value{
An object of class \code{"ar"}, see \code{\link{ar}} for details.
}
\references{

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.

Rousseeuw, P., Croux, C., Todorov, V., Ruckstuhl, A., Salibian-Barrera, M., Verbeke, T., Koller, M. and Maechler, M. (2014): robustbase: Basic Robust Statistics. R package version 0.91-1. URL \url{http://CRAN.R-project.org/package=robustbase}.
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
arrob(tss)
arrob(tss, method = "durbin-levinson")
arrob(tss, method = "robustregression")
}

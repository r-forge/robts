\encoding{utf8}
\name{arrob}
\alias{arrob}

\title{
Robust Autoregression Model
}
\description{
	Robust fit of an autogression model to a time series.  Its order is optimized with respect to the AIC or given.
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
	\item{acf.fun}{character string naming the function to calculate the autocorrelation function.  See \code{\link{acfrob}} for details.  Do not use \code{"acfrobfil"}.}
	\item{aicpenalty}{function of the model order defining a penalty term to be added to the AIC.}
}

\details{

There are many possibilities to robustly estimate an AR model

For default \code{method == "yule-walker"}, the autocorrelation function is calculated by \code{\link{acfrob}}.  Afterwards, the regression coefficients are obtained from it by solving the Yule-Walker equations, see Maronna et al. (2006), chapter 8.2 for details.

For \code{method == "durbin-levinson"}, the Durbin-Levinson algorithm is used, see Maronna et al. (2006), chapter 8.2.

For \code{method == "ols"}, the regression coefficients are computed by robustly fitting a linear model, see Maronna et al. (2006), chapter 8.2 for details.  The model fit is performed by \code{\link{lmrob}}, see Rousseeuw et al. (2014).

For \code{method == "filter"}, the regression coefficients are computed via robust filtering.  See \enc{Dürre}{Duerre} et al. (2014) and Maronna et al. (2006), chapter 8.6 for details.

For \code{method == "gm"}, generalized M estimates are used for the AR fit, see Maronna et al. (2006) for details.

See also Maronna et al. (2006), chapter 8.6.6 for optimization of the AR order.
}
\value{
An object of class \code{"ar"}, see \code{\link{ar}} for details.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

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

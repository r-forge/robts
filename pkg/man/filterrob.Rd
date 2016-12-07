\encoding{utf8}
\name{filterrob}
\alias{filterrob}
\alias{filterrob.given}

\title{
Robust Filtering of Time Series Assuming an Autoregressive Model
}
\description{
Robustly filters a given time series assuming an Autoregressive Model. The filtered series can be applied to improve estimations. Optionally, the coefficients of a fitted AR model can be used.  See chapter 8.6 of Maronna et al. (2006) and \enc{Dürre}{Duerre} et al. (2014) for details.
}
\usage{
filterrob(x, ar = NULL, method = c("fit", "given"), psifn, locfn, scalefn,
          na.action = na.fail, ...)
filterrob.given(x, ar, psifn, locfn, scalefn, na.action = na.fail)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{ar}{numeric vector of AR coefficients.}
	\item{method}{character string whether the filtering is performed using the AR coefficients given in argument \code{ar} (\code{method = "given"}) or if the AR coefficients are fitted previously by the function \code{arrob.filter} (\code{method = "fit"}).}
	\item{na.action}{function to be called to handle missing values.}
  \item{psifn}{psi function for filtering of residuals if \code{method = "pgram"}, see Details.}
	\item{locfn}{function which calculates a location estimator. Its argument must be a vector of the data. One could for example use the function \code{\link{median}}.}
  \item{scalefn}{function which calculates a scale estimator. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}}, \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
	\item{...}{additional arguments passed to function \code{arrob.filter}. This is only applicable if \code{method = "fit"}.}
}

\details{
This function filters a time series in respect of outliers. It is assumed that the time series follows an AR process with Gaussian innovations.

If the model coefficients are known or have been previously estimated, one can set \code{method = "given"} and provide the vector of AR coefficients by the argument \code{ar}. Note that this argument has the value \code{NULL} by default, which corresponds to an AR model of order zero (assuming independent innovations).

If the coefficients are not available, they can be estimated by setting \code{method = "fit"}. Estimation is done by the function \code{arrob.filter} for which further arguments can be provided by the \code{...} argument. In this case no additional information on the model fit is returned, e.g. the estimated parameters. It is difficult to assess the quality of the fit. For this reason it is recommended to fit the model in a separate step before calling this function and then use \code{method = "given"}..

Based on the model, residuals are calculated and smoothed by a psi function which is given by the argument \code{psifunc}. By default this is \code{M_psi(x, type="smooth")}, a smooth function which is the identity for absolute values smaller than 2 and zero for absolute values larger than 3. See \enc{Dürre}{Duerre} et al. (2014) for a formal definition of this psi function. 

Note that this function can be also used to handle missing values. If one replaces missing values with unusually large ones, the algorithm will replace them by best linear predictions. 
}

\value{
	Named list of the following elements:
	\item{filtered}{numeric vector with the robustly filtered version of \code{x}.}
	\item{residuals}{residuals of the fitted AR model.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): \emph{Robust Statistics: Theory and Methods}, Wiley, chapter 8.
}

\seealso{
The function \code{\link{arrob.filter}} estimates an AR model using the robust filter algorithm.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
filterrob(tss, ar = phis, method = "fit")
phis <- arrob(tss, aic = FALSE, order.max = 3)$ar
filterrob(tss, ar = phis, method = "given")
}

\keyword{Autoregressive models}

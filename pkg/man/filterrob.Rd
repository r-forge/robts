\encoding{utf8}
\name{arfilterrob}
\alias{arfilterrob}
\alias{arfilterrob.given}

\title{
Robustly Filtered Time Series
}
\description{
Filters a given time series robustly. The filtered series can be applied to improve estimations. Optionally, the coefficients of a fitted AR model can be used.  See chapter 8.6 of Maronna et al. (2006) and \enc{Dürre}{Duerre} et al. (2014) for details.
}
\usage{
filterrob(x, ar = NULL, method = c("fit", "given"), psifn, locfn, scalefn, na.action = na.fail, ...)
filterrob.given(x, ar, psifn, locfn, scalefn, na.action = na.fail)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{ar}{numeric vector AR coefficients if \code{method == "ar"}.}
	\item{method}{character string whether the filtering is performed using the AR coefficients given in argument \code{ar} or if the AR coefficients are fitted previously by the function \code{arrob.filter}.}
	\item{na.action}{function to be called to handle missing values.}
	\item{psifunc}{XXX}
	\item{locfunc}{XXX}
	\item{scalefunc}{XXX}
	\item{...}{additional arguments passed to function \code{arrob.filter}. This is only appliccable if \code{method = "fit"}.}
}

\details{
This function filters a time series in respect of outliers. It is assumed that the time series follows an AR process with Gaussian innovations.

If the model coefficients are known or have been previously estimated, one can set \code{method = "given"} and provide the vector of AR coefficients by the argument \code{ar}. Note that this argument has the value \code{NULL} by default, which corresponds to an AR model of order zero (assuming independent innovations).

If the coefficients are not available, they can be estimated by setting \code{method = "fit"}. Estimation is done by the function \code{arrob.filter} for which further arguments can be provided by the \code{...} argument. In this case no additional information on the model fit is returned, e.g. the estimated parameters. It is difficult to assess the quality of the fit. For this reason it is reccommended to fit the model in a separate step before calling this funtion und then use \code{method = "given"}..

Based on the model, residuals are calculated and smoothed by a \eqn{\Psi} function which is given by the argument \code{psifunc}. By default this is \code{\link{M_psi}(x, type="smooth")}, a smooth function which is the identity for absolute values smaller than 2 and zero for absolute values larger than 3. See \enc{Dürre}{Duerre} et al. (2014) for a formal definition of this \eqn{\Psi} function. 

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
phis <- arrob(tss, aic = FALSE, order.max = 3)$ar
filterrob(tss, ar = phis, method = "ar")
}

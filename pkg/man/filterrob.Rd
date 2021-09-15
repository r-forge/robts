\encoding{utf8}
\name{filterrob}
\alias{filterrob}
\alias{filterrob.statespace}
\alias{filterrob.recursive}
\alias{filterrob.statespaceARMA}

\title{
Robust Filtering of Time Series Assuming an Autoregressive or ARMA Model
}
\description{
Robustly filters a given time series assuming an Autoregressive or ARMA Model. The filtered series can be applied to improve estimations. Optionally, the coefficients of a fitted AR model can be used.  See Chapter 8.6 of Maronna et al. (2006) and \enc{Dürre}{Duerre} et al. (2015) for details.
}
\usage{
filterrob(x, ar, ma, var.pred, method = c("statespace", "recursive", "statespaceARMA"),
          locfn, psi.l = 2, psi.0 = 3, na.action = na.fail, ...)
filterrob.statespace(x, ar, var.pred, locfn, psi.l = 2, psi.0 = 3,
                     na.action = na.fail)
filterrob.recursive(x, ar, var.pred, locfn, psi.l = 2, psi.0 = 3,
                    na.action = na.fail)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{ar}{numeric vector of AR coefficients. If missing, the AR coefficients are estimated by robust filtering with \code{\link{arrob.filter}}.}
    \item{ma}{numeric vector of MA coefficients.}
  \item{var.pred}{numeric value with the innovation variance of the AR model. This argument needs to be provided unless argument \code{ar} is missing as well. In that case its value is taken directly from the fitted AR model.}
	\item{method}{character string specifying the method applied for robust filtering. If \code{method = "statespace"} (the default), the procedure by Maronna et al. (2006) is used. If \code{method = "recursive"}, a somewhat simpler procedure is used, see Details. If \code{method = "statespaceARMA"}, the procedure by Maronna et al. (2006) for an ARMA model is used.}
	\item{locfn}{function which calculates a location estimator for previously centering the time series. Its argument must be a vector of the data. One could for example use the function \code{\link{median}}. Ba default the mean estimation of the function \code{\link{scaleTau2}} is used.}
	\item{psi.l}{numeric value determining the psi function used for filtering, see \code{link{smoothpsi}}. For absolute values lower than \code{psi.l} the psi function is linear.}
	\item{psi.0}{numeric value determining the psi function used for filtering, see \code{link{smoothpsi}}. For absolute values larger than \code{psi.0} the psi function is zero.}
	\item{na.action}{function to be called to handle missing values.}
	\item{...}{additional arguments passed to function \code{arrob.filter}. This is only applicable if argument \code{ar} is missing.}
}

\details{
This function filters a time series in order to detect and shrink potential outliers. This filtering is based on the assumption that the time series follows an AR or ARMA process with Gaussian innovations.

If the model coefficients are known or have been previously estimated, one can provide the vector of AR coefficients by the argument \code{ar}, repectively the vector of MA coefficients by the argument \code{ma}. In that case one additionally need to specify the innovation variance by the argument \code{var.pred}. An AR model of order zero (assuming independent innovations) could be specified by \code{ar = NULL}.

If the coefficients are not provided, then in a first step an AR model is fitted to the data. Estimation is done by the function \code{arrob.filter} for which further arguments can be provided by the \code{...} argument. In this case no additional information on the model fit is returned, e.g. the estimated parameters and the innovation variance. It is difficult to assess the quality of this fit. For this reason it is recommended to fit the model by the function \code{\link{arrob}} in a separate step before calling this function and then providing the estimated coefficients and the estimated innovation variance, see Examples.

Based on the model, residuals are calculated and smoothed by a psi function which is given by the argument \code{psifunc}. By default this is \code{M_psi(x, type="smooth")}, a smooth function which is the identity for absolute values smaller than 2 and zero for absolute values larger than 3. See \code{link{smoothpsi}} and \enc{Dürre}{Duerre} et al. (2015) for a formal definition of this psi function. 

The procedure used for robust filtering can be specified by the argument \code{method} and can be one of the following possibilities:
\describe{
\item{\code{"statespace"}}{The approach proposed in Chapter 8.6 of Maronna et al. (2006). It employs a state space representation of the model.}
\item{\code{"recursive"}}{A simple strictly recursive approach, where filtering of the observations is only determined by previous but not by successive observations.}
\item{\code{"statespaceARMA"}}{The approach proposed in Chapter 8.8 of Maronna et al. (2006). It employs a state space representation of the model.}
}

Note that this function can be also used to handle missing values. If one replaces missing values with unusually large ones, the algorithm will replace them by best linear predictions. 
}

\value{
	Named list of the following elements:
	\item{filtered}{numeric vector with the robustly filtered version of \code{x}.}
	\item{residuals}{residuals of the fitted AR model.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222, \doi{10.1002/wics.1351}.

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): \emph{Robust Statistics: Theory and Methods}, Wiley, chapter 8, \doi{10.1002/0470010940}.
}

\seealso{
The function \code{\link{arrob.filter}} estimates an AR model using the robust filter algorithm.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.6), n = 100)
tss[50] <- 100 # extreme outlier in the data

# With previous model fit:
fit <- arrob(tss, method = "filter")
filterrob(tss, ar = fit$ar, var.pred = fit$var.pred, method = "recursive")

# Without previous model fit:
filterrob(tss, method = "statespace")
}

\keyword{Autoregressive models}

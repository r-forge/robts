\encoding{utf8}
\name{ARfilter}
\alias{ARfilter}

\title{
Robust Fitting of Autoregressive Models by Robust Filtering
}

\description{
Robustly fits AR models using a robust filter algorithm. See \enc{Dürre}{Duerre} et al. (2015) for details.

This is a workhorse function which is called by several other functions in this package.
}

\usage{
ARfilter(x, order.max, aicpenalty = function(p) {2*p}, psi.l = 2, psi.0 = 3)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{order.max}{integer value determining the maximum order of the AR fit.}
	\item{psi.l}{numeric value determining the used psi function. See Details.}
	\item{psi.0}{numeric value determining the used psi function. See Details.}
	\item{aicpenalty}{function of the model order, indicating the penalty for a larger model. See Details.}
	}

\details{
The function fits AR models of increasing order by a robust version of the Durbin Levinson algorithm as described in Chapter 8 of Maronna et al. (2006). The AR coefficients are estimated by minimizing a robust scale (\code{\link{scaleTau2}}) of one step ahead residuals of robustly filtered predictors. The filter process basically compares the predicted value under the estimated model with the observed one and transforms the corresponding residual with the psi function. The transformed value is then set to the sum of the predicted value and the transformed residual. For more details see Maronna et al. (2006) and \enc{Dürre}{Duerre} et al. (2015).

Following Maronna et al. (2006), the psi function should fulfill two properties. It should be the identity for small absolute values and 0 for large absolute values. Both thresholds can be determined by \code{psi.l} and \code{psi.0}. Here the psi function is chosen to be two times continuous differentiable, see \enc{Dürre}{Duerre} et al. (2015) for a formal definition.
}

\value{
List containing the following values:

\item{pacf}{numeric vector representing the successively estimated partial correlations of the model at lags 1,...,\code{p}.}
\item{var}{numeric vector containing the estimations of the residual scale for AR models of order 0,...,\code{p}.}
\item{acf}{numeric vector of estimated autocorrelations of a model of order \code{order.max} at lags 1,...,\code{p}.}
\item{aic}{numeric vector of AIC values of AR models of order 0,...,\code{p}.}
\item{filtered}{numeric matrix containing the robustly filtered time series of AR models of different orders. The first column corresponds to the AR(0) model, the second to the AR(1) model and so on.}
\item{ar}{list containing the estimated AR coefficients of AR models of different orders. The first list element corresponds to the AR(0) model, the second to the AR(1) model and so on.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

Maronna, R., Martin, D. and Yohai, V. (2006): \emph{Robust statistics}, Wiley, Chichester.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
\code{\link{arrob.filter}} for robust fitting of AR models, \code{\link{acfrob.filter}} for robust estimation of the autocorrelation function and \code{\link{filterrob}} for robust filtering of time series.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
ARfilter(tss, order.max = 5)
}

\keyword{Autoregressive models}

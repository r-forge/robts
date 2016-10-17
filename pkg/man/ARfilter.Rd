\encoding{utf8}
\name{ARfilter}
\alias{ARfilter}

\title{
Estimation of AR Models by Robust Filtering
}

\description{
Robustly estimates AR models using a robust filter algorithm. See \enc{Dürre}{Duerre} et al. (2015) for details.

This is a workhorse function which is called by several other functions in this package.
}

\usage{
ARfilter(timeseries, p, aicpenalty = function(p) {2*p}, psi.l=2, psi.0=3)
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{p}{numeric determining the order of the AR fit.}
	\item{psi.l}{numeric determining the used psi function, see Details.}
	\item{psi.0}{numeric determining the used psi function, see Details.}
	\item{aicpenalty}{function of p, indicating the penalty for a larger model, see Details.}
	}

\details{
The function fits AR models of increasing order by a robust version of the Durbin Levinson algorithm as described in Chapter 8 of Maronna et al. (2006). The AR coefficients are estimated by minimizing a robust scale (\code{\link{scaleTau2}}) of one step ahead residuals of robustly filtered predictors. The filter process basically compares the predicted value under the estimated model with the observed one and transformes the corresponding residual with the psi function. The transformed value is then set to the sum of the predicted value and the transformed residual. For more details see Maronna et al. (2006) and \enc{Dürre}{Duerre} et al. (2015).

Following Maronna et al. (2006), the psi function should fulfill two properties. It should be the identity for small absolute values and 0 for large absolute values. Both thresholds can be determined by \code{psi.l} and \code{psi.0}. Here the psi function is choosen to be two times continuous differentiable, see \enc{Dürre}{Duerre} et al. (2015) for a formal definition.
}

\value{
List containing the following values:

\item{partial autocorrelations}{numeric matrix representing the successive estimated partial correlations of the model. In the first column are the partial correlations of the AR model of order 1, in the second one the partial correlations of order 2 and so on.}
\item{variance of innovations}{numeric vector containing the estimations of the residual scale for AR models of order 1,...,\code{p}.}
\item{autocorrelation}{numeric vector of estimated autocorrelations at lags 1,...,\code{p}.}
\item{aic}{numeric vector of aic values of AR models of order 1,...,\code{p}.}
\item{robustly filtered time series}{numeric matrix containing the robustly filtered time series of AR models of different orders. The first column corresponds to the AR(1) model and so on.}
\item{AR-coefficients}{numeric matrix containing the estimated AR coefficients of AR models of different orders. The first column corresponds to the AR(1) model and the first row to the first AR coefficient and so on.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

Maronna, R., Martin, D. and Yohai, V. (2006): \emph{Robust statistics}, Wiley, Chichester.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
\code{\link{arrob}} for robust fitting of AR models, \code{\link{acfrob.filter}} for robust estimation of the autocorrelation function and \code{\link{filterrob}} for robust filtering of time series.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
ARfilter(tss, p = 5)
}
